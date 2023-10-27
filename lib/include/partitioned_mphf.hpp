#pragma once

#include "ef_sequence.hpp"
#include "quartet_wtree.hpp"
#include "mm_quartet.hpp"
#include "mphf.hpp"

namespace lphash {
namespace mphf {

class partitioned : public interface {
public:
    partitioned();
    void build(configuration const& config, std::ostream& res_strm);
    uint64_t num_bits() const noexcept;
    void print_statistics() const noexcept;

    template <typename MinimizerHasher = pthash::murmurhash2_64>
    std::vector<uint64_t> operator()(const char* contig, std::size_t length, bool streaming = true) const;

    template <typename MinimizerHasher = pthash::murmurhash2_64>
    std::vector<uint64_t> operator()(std::string const& contig, bool streaming = true) const;

    template <typename Visitor>
    void visit(Visitor& visitor);

    friend std::ostream& operator<<(std::ostream& os, const partitioned& obj);

private:
    uint64_t n_maximal, right_coll_sizes_start, none_sizes_start, none_pos_start;  // Left positions | right + coll sizes | none sizes | none positions
    quartet_wtree wtree;
    ef_sequence sizes_and_positions;

    struct mm_context_t {
        uint64_t hval;
        uint64_t global_rank;
        uint64_t local_rank;
        uint8_t type;
    };
    mm_context_t query(kmer_t kmer, uint64_t minimizer, uint32_t position) const;
    void build_inverted_index(external_memory_vector<mm_triplet_t>::const_iterator& mm_itr,
                              std::size_t number_of_distinct_minimizers);
};

template <typename MinimizerHasher>
std::vector<uint64_t> partitioned::operator()(const char* contig, std::size_t length,
                                       bool streaming) const {
    using namespace minimizer;
    std::vector<uint64_t> res;
    if (streaming) {
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
                    km[1] = (km[1] >> 2) |
                            ((static_cast<kmer_t>(3) ^ static_cast<kmer_t>(c)) << km_shift);
                    // if (canonical_m_mers && mm[0] != mm[1]) z = mm[0] < mm[1] ? 0 : 1;  //
                    // strand, if symmetric k-mer then use previous strand
                    ++nbases_since_last_break;
                    if (nbases_since_last_break >= m) {
                        current.itself = mm[z];
                        current.hash = MinimizerHasher::hash(mm[z], mm_seed).first();
                        if (buf_pos == min_pos)
                            find_brand_new_min = 1;  // old minimum out of the window
                        buffer[buf_pos] = current;
                        if (nbases_since_last_break == k) {  // first window
                            assert(buf_pos == static_cast<std::size_t>(k - m));
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
                                        mm_ctx.local_rank =
                                            fallback_kmer_order(km[z]);  // collision
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
                                for (std::size_t j = (buf_pos + 2) % buffer.size();
                                     j < buffer.size(); ++j) {
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
    } else {
        for (std::size_t i = 0; i < length - k + 1; ++i) {
            auto kmer = debug::string_to_integer_no_reverse(&contig[i], k);
            debug::triplet_t triplet =
                debug::compute_minimizer_triplet<MinimizerHasher>(kmer, k, m, mm_seed);
            uint64_t mm = triplet.first;
            uint64_t p = triplet.third;
            auto ctx = query(kmer, mm, p);
            res.push_back(ctx.hval);
        }
    }
    return res;
}

template <typename MinimizerHasher>
std::vector<uint64_t> partitioned::operator()(std::string const& contig, bool streaming) const {
    return operator()<MinimizerHasher>(contig.c_str(), contig.length(), streaming);
}

template <typename Visitor>
void partitioned::visit(Visitor& visitor) {
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

} // namespace mphf
} // namespace lphash
