#pragma once

#include "../include/constants.hpp"
#include "ef_sequence.hpp"
#include "quartet_wtree.hpp"
#include "mm_context.hpp"

namespace lphash {

template <typename MPHFType>
bool check_collisions(MPHFType const& hf, std::string const& contig, bool canonical, pthash::bit_vector_builder& population) 
{  // Note fast and dumb hashes are compared in check_streaming_correctness
    auto hashes = hf(contig, canonical, false);
    for (auto hash : hashes) {
        if (hash > hf.get_kmer_count()) {
            std::cerr << "[Error] overflow : " << hash << " > " << hf.get_kmer_count() << std::endl;
            return false;
        } else if (population.get(hash) == 1) {
            std::cerr << "[Error] collision at position (hash) : " << hash << std::endl;
            return false;
        } else
            population.set(hash);
    }
    return true;
}

template <typename MPHFType>
bool check_perfection(MPHFType const& hf, pthash::bit_vector_builder& population) 
{
    bool perfect = true;
    for (std::size_t i = 0; i < hf.get_kmer_count(); ++i)
        if (!population.get(i)) { perfect = false; }
    if (!perfect) {
        std::cerr << "[Error] Not all k-mers have been marked by a hash" << std::endl;
        return false;
    } else
        std::cerr << "[Info] Everything is ok\n";
    return perfect;
}

template <typename MPHFType>
bool check_streaming_correctness(MPHFType const& hf, std::string const& contig, bool canonical) 
{
    auto dumb_hashes = hf(contig, canonical, false);
    auto fast_hashes = hf(contig, canonical);
    if (dumb_hashes.size() != fast_hashes.size()) {
        std::cerr << "[Error] different number of hashes, maybe there were some Ns in the input "
                     "(not supported as of now)\n";
        return false;
    }
    for (std::size_t i = 0; i < dumb_hashes.size(); ++i) {
        if (dumb_hashes[i] != fast_hashes[i]) {
            std::cerr << "[Error] different hashes, maybe there were some Ns in the input (not "
                         "supported as of now)\n";
            return false;
        }
    }
    return true;
}

class mphf 
{
public:
    mphf();
    mphf(uint8_t klen, uint8_t mm_size, uint64_t seed, uint64_t total_number_of_kmers,
         uint8_t nthreads, std::string temporary_directory = "", bool verbose = false);
    void build_minimizers_mphf(std::vector<mm_triplet_t>& minimizers);
    void build_fallback_mphf(std::vector<kmer_t>& colliding_kmers);
    std::vector<uint64_t> build_inverted_index(std::vector<mm_triplet_t>& minimizers);
    uint64_t get_minimizer_L0() const noexcept;
    uint64_t get_kmer_count() const noexcept;
    uint64_t num_bits() const noexcept;

    template <typename MinimizerHasher = hash64>
    std::vector<uint64_t> operator()(const char* contig, std::size_t length,
                                     bool canonical = false) const;
    template <typename MinimizerHasher = hash64>
    std::vector<uint64_t> operator()(const char* contig, std::size_t length, bool canonical,
                                     bool dummy) const;
    template <typename MinimizerHasher = hash64>
    std::vector<uint64_t> operator()(std::string const& contig, bool canonical = false) const;
    template <typename MinimizerHasher = hash64>
    std::vector<uint64_t> operator()(std::string const& contig, bool canonical, bool dummy) const;

    template <typename MinimizerHasher = hash64>
    uint64_t barebone_streaming_query(const char* contig, std::size_t length,
                                      bool canonical = false) const;
    template <typename MinimizerHasher = hash64>
    uint64_t barebone_dumb_query(const char* contig, std::size_t length, bool canonical) const;

    void print_statistics() const noexcept;

    template <typename Visitor>
    void visit(Visitor& visitor);

    friend std::ostream& operator<<(std::ostream& os, const mphf& obj);

private:
    pthash::build_configuration mphf_configuration;
    uint8_t k, m;
    uint64_t mm_seed;
    uint64_t nkmers;
    uint64_t distinct_minimizers;
    uint64_t n_maximal, right_coll_sizes_start, none_sizes_start, none_pos_start;  // Left positions | right + coll sizes | none sizes | none positions
    pthash_mphf_t minimizer_order;
    pthash_mphf_t fallback_kmer_order;
    quartet_wtree wtree;
    ef_sequence<true> sizes_and_positions;

    struct mm_context_t {
        uint64_t hval;
        uint64_t global_rank;
        uint64_t local_rank;
        uint8_t type;
    };
    mm_context_t query(kmer_t kmer, uint64_t minimizer, uint32_t position) const;
};

template <typename MinimizerHasher>
std::vector<uint64_t> mphf::operator()(const char* contig, std::size_t length, [[maybe_unused]] bool canonical_m_mers) const 
{
    using namespace minimizer;
    std::vector<uint64_t> res;
    if (length < k) return res;
    res.reserve(length - k + 1);
    uint64_t shift = 2 * (m - 1);
    uint64_t mask = (1ULL << (2 * m)) - 1;
    uint64_t mm[2] = {0, 0};
    uint64_t km_shift = 2 * (k - 1);
    kmer_t km_mask = (static_cast<kmer_t>(1) << (2 * k)) - 1;
    kmer_t km[2] = {0, 0};
    uint64_t nbases_since_last_break = 0;
    uint8_t find_brand_new_min = 0;
    uint32_t p1 = 0;

    std::vector<mm_quartet_t> buffer(k - m + 1);
    std::size_t buf_pos, min_pos;
    mm_quartet_t current;
    mm_context_t mm_ctx;
    int c;
    uint8_t z;

    assert(k >= m);

    buf_pos = 0;
    min_pos = buffer.size();
    z = 0;
    for (uint64_t i = 0; i < length; ++i) {
        c = constants::seq_nt4_table[static_cast<uint8_t>(contig[i])];
        if (c < 4) [[likely]] {
                mm[0] = (mm[0] << 2 | c) & mask;            /* forward m-mer */
                mm[1] = (mm[1] >> 2) | (3ULL ^ c) << shift; /* reverse m-mer */
                km[0] = (km[0] << 2 | static_cast<kmer_t>(c)) & km_mask;
                km[1] =
                    (km[1] >> 2) | ((static_cast<kmer_t>(3) ^ static_cast<kmer_t>(c)) << km_shift);
                if (canonical_m_mers && mm[0] != mm[1])
                    z = mm[0] < mm[1] ? 0 : 1;  // strand, if symmetric k-mer then use previous strand
                ++nbases_since_last_break;
                if (nbases_since_last_break >= m) {
                    current.itself = mm[z];
                    current.hash = MinimizerHasher::hash(mm[z], mm_seed).first();
                    if (buf_pos == min_pos)
                        find_brand_new_min = 1;  // old minimum out of the window
                    buffer[buf_pos] = current;
                    if (nbases_since_last_break == k) {  // first window
                        assert(buf_pos == k - m);
                        find_brand_new_min = 1;
                    } else if (nbases_since_last_break > k) {
                        if (buffer[min_pos].hash > buffer[buf_pos].hash) {  // new minimum
                            p1 = k - m;
                            min_pos = buf_pos;
                            find_brand_new_min = 2;
                        }
                    }
                    switch (find_brand_new_min) {
                        case 0: {
                            if (nbases_since_last_break >= k) {
                                if (mm_ctx.type == (NONE + 1)) {
                                    mm_ctx.local_rank = fallback_kmer_order(km[z]);  // collision
                                } else if (mm_ctx.type == RIGHT_OR_COLLISION ||
                                           mm_ctx.type == NONE) {
                                    ++mm_ctx.local_rank;
                                } else {
                                    --mm_ctx.local_rank;
                                }
                                mm_ctx.hval = mm_ctx.global_rank + mm_ctx.local_rank;
                                res.push_back(mm_ctx.hval);
                            }
                        } break;
                        case 1: {
                            min_pos = (buf_pos + 1) % buffer.size();
                            p1 = 0;
                            uint32_t tmp = 1;
                            for (std::size_t j = (buf_pos + 2) % buffer.size(); j < buffer.size();
                                 ++j) {
                                if (buffer[min_pos].hash > buffer[j].hash) {
                                    min_pos = j;
                                    p1 = tmp;
                                }
                                ++tmp;
                            }
                            for (std::size_t j = 0; j <= (buf_pos + 2) % buffer.size(); ++j) {
                                if (buffer[min_pos].hash > buffer[j].hash) {
                                    min_pos = j;
                                    p1 = tmp;
                                }
                                ++tmp;
                            }
                        }
                            [[fallthrough]];
                        case 2:
                            find_brand_new_min = 0;
                            mm_ctx = query(km[z], buffer[min_pos].itself, p1);
                            res.push_back(mm_ctx.hval);
                            break;
                        default:
                            throw std::runtime_error(
                                "[Error] mphf operator(): this should never happen");
                    }
                    buf_pos = (buf_pos + 1) % buffer.size();
                }
            }
        else {
            // std::cerr << "reset\n";
            nbases_since_last_break = 0;
            buf_pos = 0;
        }
    }
    return res;
}

template <typename MinimizerHasher>
std::vector<uint64_t> mphf::operator()(const char* contig, std::size_t length,
                                       [[maybe_unused]] bool canonical_m_mers,
                                       [[maybe_unused]] bool dummy) const {
    std::vector<uint64_t> res;
    for (std::size_t i = 0; i < length - k + 1; ++i) {
        auto kmer = debug::string_to_integer_no_reverse(&contig[i], k);
        debug::triplet_t triplet =
            debug::compute_minimizer_triplet<MinimizerHasher>(kmer, k, m, mm_seed);
        uint64_t mm = triplet.first;
        uint64_t p = triplet.third;
        auto ctx = query(kmer, mm, p);
        res.push_back(ctx.hval);
    }
    return res;
}

template <typename MinimizerHasher>
std::vector<uint64_t> mphf::operator()(std::string const& contig,
                                       [[maybe_unused]] bool canonical_m_mers) const {
    return operator()<MinimizerHasher>(contig.c_str(), contig.length(), canonical_m_mers);
}

template <typename MinimizerHasher>
std::vector<uint64_t> mphf::operator()(std::string const& contig,
                                       [[maybe_unused]] bool canonical_m_mers,
                                       [[maybe_unused]] bool dummy) const {
    return operator()<MinimizerHasher>(contig.c_str(), contig.length(), canonical_m_mers, dummy);
}

template <typename MinimizerHasher>
uint64_t mphf::barebone_streaming_query(const char* contig, std::size_t length,
                                        [[maybe_unused]] bool canonical_m_mers) const {
    using namespace minimizer;
    uint64_t res = 0;
    if (length < k) return res;
    uint64_t shift = 2 * (m - 1);
    uint64_t mask = (1ULL << (2 * m)) - 1;
    uint64_t mm[2] = {0, 0};
    uint64_t km_shift = 2 * (k - 1);
    kmer_t km_mask = (static_cast<kmer_t>(1) << (2 * k)) - 1;
    kmer_t km[2] = {0, 0};
    uint64_t nbases_since_last_break = 0;
    uint8_t find_brand_new_min = 0;
    uint32_t p1 = 0;

    std::vector<mm_quartet_t> buffer(k - m + 1);
    std::size_t buf_pos, min_pos;
    mm_quartet_t current;
    mm_context_t mm_ctx;
    int c;
    uint8_t z;

    assert(k >= m);

    buf_pos = 0;
    min_pos = buffer.size();
    z = 0;
    for (uint64_t i = 0; i < length; ++i) {
        c = constants::seq_nt4_table[static_cast<uint8_t>(contig[i])];
        if (c < 4) [[likely]] {
                mm[0] = (mm[0] << 2 | c) & mask;            /* forward m-mer */
                mm[1] = (mm[1] >> 2) | (3ULL ^ c) << shift; /* reverse m-mer */
                km[0] = (km[0] << 2 | static_cast<kmer_t>(c)) & km_mask;
                km[1] =
                    (km[1] >> 2) | ((static_cast<kmer_t>(3) ^ static_cast<kmer_t>(c)) << km_shift);
                if (canonical_m_mers && mm[0] != mm[1])
                    z = mm[0] < mm[1] ? 0
                                      : 1;  // strand, if symmetric k-mer then use previous strand
                ++nbases_since_last_break;
                if (nbases_since_last_break >= m) {
                    current.itself = mm[z];
                    current.hash = MinimizerHasher::hash(mm[z], mm_seed).first();
                    if (buf_pos == min_pos)
                        find_brand_new_min = 1;  // old minimum out of the window
                    buffer[buf_pos] = current;
                    if (nbases_since_last_break == k) {  // first window
                        assert(buf_pos == k - m);
                        find_brand_new_min = 1;
                    } else if (nbases_since_last_break > k) {
                        if (buffer[min_pos].hash > buffer[buf_pos].hash) {  // new minimum
                            p1 = k - m;
                            min_pos = buf_pos;
                            find_brand_new_min = 2;
                        }
                    }
                    switch (find_brand_new_min) {
                        case 0: {
                            if (nbases_since_last_break >= k) {
                                if (mm_ctx.type == (NONE + 1)) {
                                    mm_ctx.local_rank = fallback_kmer_order(km[z]);  // collision
                                } else if (mm_ctx.type == RIGHT_OR_COLLISION ||
                                           mm_ctx.type == NONE) {
                                    ++mm_ctx.local_rank;
                                } else {
                                    --mm_ctx.local_rank;
                                }
                                mm_ctx.hval = mm_ctx.global_rank + mm_ctx.local_rank;
                                ++res;
                            }
                        } break;
                        case 1: {
                            min_pos = (buf_pos + 1) % buffer.size();
                            p1 = 0;
                            uint32_t tmp = 1;
                            for (std::size_t j = (buf_pos + 2) % buffer.size(); j < buffer.size();
                                 ++j) {
                                if (buffer[min_pos].hash > buffer[j].hash) {
                                    min_pos = j;
                                    p1 = tmp;
                                }
                                ++tmp;
                            }
                            for (std::size_t j = 0; j <= (buf_pos + 2) % buffer.size(); ++j) {
                                if (buffer[min_pos].hash > buffer[j].hash) {
                                    min_pos = j;
                                    p1 = tmp;
                                }
                                ++tmp;
                            }
                        }
                            [[fallthrough]];
                        case 2:
                            find_brand_new_min = 0;
                            mm_ctx = query(km[z], buffer[min_pos].itself, p1);
                            ++res;
                            break;
                        default:
                            throw std::runtime_error(
                                "[Error] mphf operator(): this should never happen");
                    }
                    buf_pos = (buf_pos + 1) % buffer.size();
                }
            }
        else {
            nbases_since_last_break = 0;
            buf_pos = 0;
        }
    }
    return res;
}

template <typename MinimizerHasher>
uint64_t mphf::barebone_dumb_query(const char* contig, std::size_t length,
                                   [[maybe_unused]] bool canonical_m_mers) const {
    uint64_t res = 0;
    for (std::size_t i = 0; i < length - k + 1; ++i) {
        auto kmer = debug::string_to_integer_no_reverse(&contig[i], k);
        debug::triplet_t triplet =
            debug::compute_minimizer_triplet<MinimizerHasher>(kmer, k, m, mm_seed);
        uint64_t mm = triplet.first;
        uint64_t p = triplet.third;
        [[maybe_unused]] auto ctx = query(kmer, mm, p);
        // essentials::do_not_optimize_away(ctx);
        ++res;
    }
    return res;
}

template <typename Visitor>
void mphf::visit(Visitor& visitor) {
    visitor.visit(k);
    visitor.visit(m);
    visitor.visit(mm_seed);
    visitor.visit(nkmers);
    visitor.visit(distinct_minimizers);
    visitor.visit(n_maximal);
    visitor.visit(right_coll_sizes_start);
    visitor.visit(none_sizes_start);
    visitor.visit(none_pos_start);
    visitor.visit(minimizer_order);
    visitor.visit(wtree);
    visitor.visit(sizes_and_positions);
    visitor.visit(fallback_kmer_order);
}

// bool check_collisions(mphf const& hf, std::string const& contig, bool canonical, pthash::bit_vector_builder& population); 
// bool check_perfection(mphf const& hf, pthash::bit_vector_builder& population); 
// bool check_streaming_correctness(mphf const& hf, std::string const& contig, bool canonical);

class mphf_alt {
public:
    mphf_alt();
    mphf_alt(uint8_t klen, uint8_t mm_size, uint64_t seed, uint64_t total_number_of_kmers,
             uint8_t nthreads, std::string temporary_directory = "", bool verbose = false);
    std::vector<uint64_t> build_index(std::vector<mm_triplet_t>& minimizers);
    void build_fallback_mphf(std::vector<kmer_t>& colliding_kmers);
    uint64_t get_minimizer_L0() const noexcept;
    uint64_t get_kmer_count() const noexcept;
    uint64_t num_bits() const noexcept;

    template <typename MinimizerHasher = hash64>
    std::vector<uint64_t> operator()(const char* contig, std::size_t length,
                                     bool canonical = false) const;
    template <typename MinimizerHasher = hash64>
    std::vector<uint64_t> operator()(const char* contig, std::size_t length, bool canonical,
                                     bool dummy) const;
    template <typename MinimizerHasher = hash64>
    std::vector<uint64_t> operator()(std::string const& contig, bool canonical = false) const;
    template <typename MinimizerHasher = hash64>
    std::vector<uint64_t> operator()(std::string const& contig, bool canonical, bool dummy) const;

    template <typename MinimizerHasher = hash64>
    uint64_t barebone_streaming_query(const char* contig, std::size_t length,
                                      bool canonical = false) const;
    template <typename MinimizerHasher = hash64>
    uint64_t barebone_dumb_query(const char* contig, std::size_t length, bool canonical) const;

    void print_statistics() const noexcept;

    template <typename Visitor>
    void visit(Visitor& visitor);

    friend std::ostream& operator<<(std::ostream& os, const mphf_alt& obj);

private:
    pthash::build_configuration mphf_configuration;
    uint8_t k, m;
    uint64_t mm_seed;
    uint64_t nkmers;
    uint64_t distinct_minimizers;
    uint64_t num_kmers_in_main_index;
    pthash_mphf_t minimizer_order;
    ef_sequence<true> positions;
    ef_sequence<true> sizes;
    pthash_mphf_t fallback_kmer_order;

    struct mm_context_t {
        uint64_t hval;
        uint64_t collision;
    };
    mm_context_t query(kmer_t kmer, uint64_t minimizer, uint32_t position) const;
};

template <typename MinimizerHasher>
std::vector<uint64_t> mphf_alt::operator()(const char* contig, std::size_t length,
                                           [[maybe_unused]] bool canonical_m_mers) const {
    using namespace minimizer;
    std::vector<uint64_t> res;
    if (length < k) return res;
    res.reserve(length - k + 1);
    uint64_t shift = 2 * (m - 1);
    uint64_t mask = (1ULL << (2 * m)) - 1;
    uint64_t mm[2] = {0, 0};
    uint64_t km_shift = 2 * (k - 1);
    kmer_t km_mask = (static_cast<kmer_t>(1) << (2 * k)) - 1;
    kmer_t km[2] = {0, 0};
    uint64_t nbases_since_last_break = 0;
    uint8_t find_brand_new_min = 0;
    uint32_t p1 = 0;

    std::vector<mm_quartet_t> buffer(k - m + 1);
    std::size_t buf_pos, min_pos;
    mm_quartet_t current;
    mm_context_t mm_ctx;
    int c;
    uint8_t z;

    assert(k >= m);

    buf_pos = 0;
    min_pos = buffer.size();
    z = 0;
    for (uint64_t i = 0; i < length; ++i) {
        c = constants::seq_nt4_table[static_cast<uint8_t>(contig[i])];
        if (c < 4) [[likely]] {
                mm[0] = (mm[0] << 2 | c) & mask;            /* forward m-mer */
                mm[1] = (mm[1] >> 2) | (3ULL ^ c) << shift; /* reverse m-mer */
                km[0] = (km[0] << 2 | static_cast<kmer_t>(c)) & km_mask;
                km[1] =
                    (km[1] >> 2) | ((static_cast<kmer_t>(3) ^ static_cast<kmer_t>(c)) << km_shift);
                if (canonical_m_mers && mm[0] != mm[1])
                    z = mm[0] < mm[1] ? 0
                                      : 1;  // strand, if symmetric k-mer then use previous strand
                ++nbases_since_last_break;
                if (nbases_since_last_break >= m) {
                    current.itself = mm[z];
                    current.hash = MinimizerHasher::hash(mm[z], mm_seed).first();
                    if (buf_pos == min_pos)
                        find_brand_new_min = 1;  // old minimum out of the window
                    buffer[buf_pos] = current;
                    if (nbases_since_last_break == k) {  // first window
                        assert(buf_pos == k - m);
                        find_brand_new_min = 1;
                    } else if (nbases_since_last_break > k) {
                        if (buffer[min_pos].hash > buffer[buf_pos].hash) {  // new minimum
                            p1 = k - m;
                            min_pos = buf_pos;
                            find_brand_new_min = 2;
                        }
                    }
                    switch (find_brand_new_min) {
                        case 0: {
                            if (nbases_since_last_break >= k) {
                                if (mm_ctx.collision)
                                    mm_ctx.hval = fallback_kmer_order(km[z]) +
                                                  sizes.access(sizes.size() - 1);  // collision
                                else
                                    ++mm_ctx.hval;
                                res.push_back(mm_ctx.hval);
                            }
                        } break;
                        case 1: {
                            min_pos = (buf_pos + 1) % buffer.size();
                            p1 = 0;
                            uint32_t tmp = 1;
                            for (std::size_t j = (buf_pos + 2) % buffer.size(); j < buffer.size();
                                 ++j) {
                                if (buffer[min_pos].hash > buffer[j].hash) {
                                    min_pos = j;
                                    p1 = tmp;
                                }
                                ++tmp;
                            }
                            for (std::size_t j = 0; j <= (buf_pos + 2) % buffer.size(); ++j) {
                                if (buffer[min_pos].hash > buffer[j].hash) {
                                    min_pos = j;
                                    p1 = tmp;
                                }
                                ++tmp;
                            }
                        }
                            [[fallthrough]];
                        case 2:
                            find_brand_new_min = 0;
                            mm_ctx = query(km[z], buffer[min_pos].itself, p1);
                            res.push_back(mm_ctx.hval);
                            break;
                        default:
                            throw std::runtime_error(
                                "[Error] mphf operator(): this should never happen");
                    }
                    buf_pos = (buf_pos + 1) % buffer.size();
                }
            }
        else {
            // std::cerr << "reset\n";
            nbases_since_last_break = 0;
            buf_pos = 0;
        }
    }
    return res;
}

template <typename MinimizerHasher>
std::vector<uint64_t> mphf_alt::operator()(const char* contig, std::size_t length,
                                           [[maybe_unused]] bool canonical_m_mers,
                                           [[maybe_unused]] bool dummy) const {
    std::vector<uint64_t> res;
    for (std::size_t i = 0; i < length - k + 1; ++i) {
        auto kmer = debug::string_to_integer_no_reverse(&contig[i], k);
        debug::triplet_t triplet =
            debug::compute_minimizer_triplet<MinimizerHasher>(kmer, k, m, mm_seed);
        uint64_t mm = triplet.first;
        uint64_t p = triplet.third;
        auto ctx = query(kmer, mm, p);
        res.push_back(ctx.hval);
        // std::cerr << ctx.hval << "\n";
    }
    return res;
}

template <typename MinimizerHasher>
std::vector<uint64_t> mphf_alt::operator()(std::string const& contig,
                                           [[maybe_unused]] bool canonical_m_mers) const {
    return operator()<MinimizerHasher>(contig.c_str(), contig.length(), canonical_m_mers);
}

template <typename MinimizerHasher>
std::vector<uint64_t> mphf_alt::operator()(std::string const& contig,
                                           [[maybe_unused]] bool canonical_m_mers,
                                           [[maybe_unused]] bool dummy) const {
    return operator()<MinimizerHasher>(contig.c_str(), contig.length(), canonical_m_mers, dummy);
}

template <typename MinimizerHasher>
uint64_t mphf_alt::barebone_streaming_query(const char* contig, std::size_t length,
                                            [[maybe_unused]] bool canonical_m_mers) const {
    using namespace minimizer;
    uint64_t res = 0;
    if (length < k) return res;
    uint64_t shift = 2 * (m - 1);
    uint64_t mask = (1ULL << (2 * m)) - 1;
    uint64_t mm[2] = {0, 0};
    uint64_t km_shift = 2 * (k - 1);
    kmer_t km_mask = (static_cast<kmer_t>(1) << (2 * k)) - 1;
    kmer_t km[2] = {0, 0};
    uint64_t nbases_since_last_break = 0;
    uint8_t find_brand_new_min = 0;
    uint32_t p1 = 0;

    std::vector<mm_quartet_t> buffer(k - m + 1);
    std::size_t buf_pos, min_pos;
    mm_quartet_t current;
    mm_context_t mm_ctx;
    int c;
    uint8_t z;

    assert(k >= m);

    buf_pos = 0;
    min_pos = buffer.size();
    z = 0;
    for (uint64_t i = 0; i < length; ++i) {
        c = constants::seq_nt4_table[static_cast<uint8_t>(contig[i])];
        if (c < 4) [[likely]] {
                mm[0] = (mm[0] << 2 | c) & mask;            /* forward m-mer */
                mm[1] = (mm[1] >> 2) | (3ULL ^ c) << shift; /* reverse m-mer */
                km[0] = (km[0] << 2 | static_cast<kmer_t>(c)) & km_mask;
                km[1] =
                    (km[1] >> 2) | ((static_cast<kmer_t>(3) ^ static_cast<kmer_t>(c)) << km_shift);
                if (canonical_m_mers && mm[0] != mm[1])
                    z = mm[0] < mm[1] ? 0
                                      : 1;  // strand, if symmetric k-mer then use previous strand
                ++nbases_since_last_break;
                if (nbases_since_last_break >= m) {
                    current.itself = mm[z];
                    current.hash = MinimizerHasher::hash(mm[z], mm_seed).first();
                    if (buf_pos == min_pos)
                        find_brand_new_min = 1;  // old minimum out of the window
                    buffer[buf_pos] = current;
                    if (nbases_since_last_break == k) {  // first window
                        assert(buf_pos == k - m);
                        find_brand_new_min = 1;
                    } else if (nbases_since_last_break > k) {
                        if (buffer[min_pos].hash > buffer[buf_pos].hash) {  // new minimum
                            p1 = k - m;
                            min_pos = buf_pos;
                            find_brand_new_min = 2;
                        }
                    }
                    switch (find_brand_new_min) {
                        case 0: {
                            if (nbases_since_last_break >= k) {
                                if (mm_ctx.collision)
                                    mm_ctx.hval = fallback_kmer_order(km[z]) +
                                                  sizes.access(sizes.size() - 1);  // collision
                                else
                                    ++mm_ctx.hval;
                                ++res;
                            }
                        } break;
                        case 1: {
                            min_pos = (buf_pos + 1) % buffer.size();
                            p1 = 0;
                            uint32_t tmp = 1;
                            for (std::size_t j = (buf_pos + 2) % buffer.size(); j < buffer.size();
                                 ++j) {
                                if (buffer[min_pos].hash > buffer[j].hash) {
                                    min_pos = j;
                                    p1 = tmp;
                                }
                                ++tmp;
                            }
                            for (std::size_t j = 0; j <= (buf_pos + 2) % buffer.size(); ++j) {
                                if (buffer[min_pos].hash > buffer[j].hash) {
                                    min_pos = j;
                                    p1 = tmp;
                                }
                                ++tmp;
                            }
                        }
                            [[fallthrough]];
                        case 2:
                            find_brand_new_min = 0;
                            mm_ctx = query(km[z], buffer[min_pos].itself, p1);
                            ++res;
                            break;
                        default:
                            throw std::runtime_error(
                                "[Error] mphf operator(): this should never happen");
                    }
                    buf_pos = (buf_pos + 1) % buffer.size();
                }
            }
        else {
            // std::cerr << "reset\n";
            nbases_since_last_break = 0;
            buf_pos = 0;
        }
    }
    return res;
}

template <typename MinimizerHasher>
uint64_t mphf_alt::barebone_dumb_query(const char* contig, std::size_t length,
                                       [[maybe_unused]] bool canonical_m_mers) const {
    uint64_t res = 0;
    for (std::size_t i = 0; i < length - k + 1; ++i) {
        auto kmer = debug::string_to_integer_no_reverse(&contig[i], k);
        debug::triplet_t triplet =
            debug::compute_minimizer_triplet<MinimizerHasher>(kmer, k, m, mm_seed);
        uint64_t mm = triplet.first;
        uint64_t p = triplet.third;
        [[maybe_unused]] auto ctx = query(kmer, mm, p);
        ++res;
    }
    return res;
}

template <typename Visitor>
void mphf_alt::visit(Visitor& visitor) {
    visitor.visit(k);
    visitor.visit(m);
    visitor.visit(mm_seed);
    visitor.visit(nkmers);
    visitor.visit(distinct_minimizers);
    visitor.visit(num_kmers_in_main_index);
    visitor.visit(minimizer_order);
    visitor.visit(positions);
    visitor.visit(sizes);
    visitor.visit(fallback_kmer_order);
}

}  // namespace lphash