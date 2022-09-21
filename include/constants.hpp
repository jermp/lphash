#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <array>
#include <string>
#include <iostream>

#include "../external/pthash/include/pthash.hpp"

#include <bitset>

namespace lphash {
namespace constants {

constexpr uint64_t max_k = 31;  // max *odd* size that can be packed into 64 bits
constexpr uint64_t invalid_uint64 = uint64_t(-1);
constexpr uint32_t invalid_uint32 = uint32_t(-1);
constexpr uint64_t seed = 1;
constexpr double c = 3.0;  // for PTHash
constexpr uint64_t min_l = 6;
constexpr uint64_t max_l = 12;
static const std::string default_tmp_dirname(".");
constexpr bool forward_orientation = 0;
constexpr bool backward_orientation = 1;

constexpr std::array<uint8_t, 256> seq_nt4_table = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

}  // namespace constants

struct murmurhash2_64 {
    // specialization for uint64_t
    static inline uint64_t hash(uint64_t val, uint64_t seed) {
        return pthash::MurmurHash2_64(reinterpret_cast<char const*>(&val), sizeof(val), seed);
    }
};

struct build_configuration {
    build_configuration()
        : k(31)
        , m(17)
        , seed(constants::seed)

        , l(constants::min_l)
        , c(constants::c)

        , canonical_parsing(false)
        , weighted(false)
        , verbose(true)

        , tmp_dirname(constants::default_tmp_dirname) {}

    uint64_t k;  // kmer size
    uint64_t m;  // minimizer size
    uint64_t seed;

    uint64_t l;  // drive dictionary trade-off
    double c;    // drive PTHash trade-off

    bool canonical_parsing;
    bool weighted;
    bool verbose;

    std::string tmp_dirname;

    void print() const {
        std::cout << "k = " << k << ", m = " << m << ", seed = " << seed << ", l = " << l
                  << ", c = " << c
                  << ", canonical_parsing = " << (canonical_parsing ? "true" : "false")
                  << ", weighted = " << (weighted ? "true" : "false") << std::endl;
    }
};

//------------------------------------------------------------------------------------------------------------------

namespace debug {

struct triplet_t {
    uint64_t first, second, third;
};

// static uint64_t char_to_uint64(char c) { return (c >> 1) & 3; }
static uint64_t char_to_uint64(char c) { return constants::seq_nt4_table[static_cast<uint8_t>(c)] & 3; }

[[maybe_unused]] static uint64_t string_to_uint64_no_reverse(char const* str, uint64_t k) {
    assert(k <= 32);
    uint64_t y;
    y = 0;
    for (uint64_t i = 0; i != k; ++i) {
        //x += char_to_uint64(str[i]) << (2 * i);
        y = (y << 2) | char_to_uint64(str[i]);
        //std::cerr << x << " " << y << "\n"; 
        //assert(x == y);
    }
    return y;
}

template <typename Hasher = lphash::murmurhash2_64>
static void print_hashes(std::string contig, uint64_t m, uint64_t seed) {
    for (uint64_t i = 0; i < contig.length() - m + 1; ++i) {
        uint64_t pmmer = string_to_uint64_no_reverse(&contig.data()[i], m);
        uint64_t hash = Hasher::hash(pmmer, seed);
        std::cerr << "[" << i << "] : " << pmmer << " " << hash << "\n";
    }
}

template <typename Hasher = lphash::murmurhash2_64>
static triplet_t compute_minimizer_triplet(uint64_t kmer, uint64_t k, uint64_t m, uint64_t seed) {
    assert(m < 32);
    assert(m <= k);
    uint64_t min_hash = uint64_t(-1);
    uint64_t minimizer = uint64_t(-1);
    uint64_t mask = (uint64_t(1) << (2 * m)) - 1;
    uint64_t pos = 0;
    for (uint64_t i = 0; i != k - m + 1; ++i) {
        uint64_t sub_kmer = kmer & mask;
        uint64_t hash = Hasher::hash(sub_kmer, seed);
        if (hash < min_hash) {
            min_hash = hash;
            minimizer = sub_kmer;
            pos = i;
        }
        kmer >>= 2;
    }
    return triplet_t{minimizer, min_hash, k-(pos+m)};
}

template <typename Hasher = lphash::murmurhash2_64>
static void compute_minimizers_naive(std::string const& contig, uint64_t k, uint64_t m, uint64_t seed) {
    triplet_t prev = {uint64_t(-1), uint64_t(-1), uint64_t(-1)};
    for (std::size_t i = 0; i < contig.size() - k + 1; ++i) {
        uint64_t uint64_kmer = string_to_uint64_no_reverse(&contig.data()[i], k);
        //std::cerr << std::bitset<64>(uint64_kmer) << "\n";
        auto curr = compute_minimizer_triplet(uint64_kmer, k, m, seed);
        if (prev.first != curr.first or (prev.third - curr.third > 1)) std::cerr << curr.first << " " << curr.second << " " << curr.third << "\n";
        prev = curr;
    }
}

} // namespace debug

} // namespace lphash

#endif