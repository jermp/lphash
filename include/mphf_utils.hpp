#pragma once

#include "constants.hpp"
#include "external_memory_vector.hpp"
#include "util.hpp"

namespace lphash {

class mm_itr_t {
public:
    mm_itr_t(external_memory_vector<mm_triplet_t, false>::const_iterator& mm_itr);
    void operator++();
    uint64_t operator*() const;

private:
    external_memory_vector<mm_triplet_t, false>::const_iterator& m_iterator;
};

class km_itr_t {
public:
    km_itr_t(external_memory_vector<kmer_t, false>::const_iterator& mm_itr);
    void operator++();
    kmer_t const& operator*() const;

private:
    external_memory_vector<kmer_t, false>::const_iterator& m_iterator;
};

class pos_itr_t {
public:
    typedef uint8_t value_type;
    pos_itr_t(external_memory_vector<mm_triplet_t>::const_iterator& mm_itr);
    void operator++();
    uint8_t const& operator*() const;

private:
    external_memory_vector<mm_triplet_t>::const_iterator& m_iterator;
};

class size_itr_t {
public:
    typedef uint8_t value_type;
    size_itr_t(external_memory_vector<mm_triplet_t>::const_iterator& mm_itr);
    void operator++();
    uint8_t const& operator*() const;

private:
    external_memory_vector<mm_triplet_t>::const_iterator& m_iterator;
};

template <typename MPHF>
bool check_collisions(
    MPHF const& f, char const* contig, std::size_t contig_len,
    pthash::bit_vector_builder&
        population) {  // Note fast and dumb hashes are compared in check_streaming_correctness
    auto hashes = f(contig, contig_len, false);
    for (auto hash : hashes) {
        if (hash > f.get_kmer_count()) {
            std::cerr << "[Error] overflow : " << hash << " > " << f.get_kmer_count() << std::endl;
            return false;
        } else if (population.get(hash) == 1) {
            std::cerr << "[Error] collision at position (hash) : " << hash << std::endl;
            return false;
        } else
            population.set(hash);
    }
    return true;
}

template <typename MPHF>
bool check_perfection(MPHF const& f, pthash::bit_vector_builder& population) {
    bool perfect = true;
    for (std::size_t i = 0; i < f.get_kmer_count(); ++i)
        if (!population.get(i)) { perfect = false; }
    if (!perfect) {
        std::cerr << "[Error] Not all k-mers have been marked by a hash" << std::endl;
        return false;
    } else
        std::cerr << "[Info] Everything is ok\n";
    return perfect;
}

template <typename MPHF>
bool check_streaming_correctness(MPHF const& f, char const* contig, std::size_t contig_len) {
    auto dumb_hashes = f(contig, contig_len, false);
    auto fast_hashes = f(contig, contig_len);
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

namespace debug {

struct triplet_t {
    uint64_t first, second, third;
};

static kmer_t char_to_uint(char c) { return constants::seq_nt4_table[static_cast<uint8_t>(c)] & 3; }

[[maybe_unused]] static kmer_t string_to_integer_no_reverse(const char* str, uint64_t k) {
    assert(k <= 64);
    kmer_t y;
    y = 0;
    for (uint64_t i = 0; i != k; ++i) { y = (y << 2) | char_to_uint(str[i]); }
    return y;
}

template <typename Hasher = pthash::murmurhash2_64>
static triplet_t compute_minimizer_triplet(kmer_t kmer, uint64_t k, uint64_t m, uint64_t seed) {
    assert(m <= 32);
    assert(m <= k);
    uint64_t min_hash = uint64_t(-1);
    kmer_t minimizer = kmer_t(-1);
    kmer_t mask = (kmer_t(1) << (2 * m)) - 1;
    uint64_t pos = 0;
    for (uint64_t i = 0; i != k - m + 1; ++i) {
        kmer_t mmer = kmer & mask;
        uint64_t hash = Hasher::hash(mmer, seed).first();
        if (hash <= min_hash) {
            min_hash = hash;
            minimizer = mmer;
            pos = i;
        }
        kmer = kmer >> 2;
    }
    return triplet_t{static_cast<uint64_t>(minimizer), min_hash, k - (pos + m)};
}

}  // namespace debug

}  // namespace lphash