#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <array>
#include <string>
#include <iostream>

#include "../external/pthash/include/pthash.hpp"

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

extern const std::array<uint8_t, 256> seq_nt4_table; // = {
//     4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//     4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//     4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//     4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//     4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//     4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//     4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
//     4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
// };

}

struct mm_triplet_t {
    uint64_t itself;
    uint32_t p1;
    uint32_t size;
};

struct kmer128_t {
    kmer128_t() {upper = 0; lower = 0;};
    kmer128_t(uint8_t v) {upper = 0; lower = static_cast<uint64_t>(v);};
    kmer128_t(int v) {upper = 0; lower = static_cast<uint64_t>(v);};
    kmer128_t(uint64_t v) {upper = 0; lower = v;};
    // operator uint64_t() const {return lower;};
    uint64_t upper;
    uint64_t lower;
};

// typedef kmer128_t kmer_t;
typedef uint64_t kmer_t;

struct hash64 : public pthash::murmurhash2_64 {
    // specialization for kmer128_t
    static inline pthash::hash64 hash(kmer128_t val, uint64_t seed) {
        return pthash::MurmurHash2_64(reinterpret_cast<char const*>(&val), sizeof(val), seed);
    }

    // generic range of bytes
    static inline pthash::hash64 hash(pthash::byte_range range, uint64_t seed) {
        return murmurhash2_64::hash(range, seed);
    }

    // specialization for std::string
    static inline pthash::hash64 hash(std::string const& val, uint64_t seed) {
        return murmurhash2_64::hash(val, seed);
    }

    // specialization for uint64_t
    static inline pthash::hash64 hash(uint64_t val, uint64_t seed) {
        return murmurhash2_64::hash(val, seed);
    }
};

typedef pthash::single_phf<hash64, pthash::dictionary_dictionary, true> pthash_mphf_t;

bool operator< (mm_triplet_t const& a, mm_triplet_t const& b);

std::ostream &operator<<(std::ostream &os, kmer128_t const& val);

bool operator< (kmer128_t const& a, kmer128_t const& b);

bool operator== (kmer128_t const& a, kmer128_t const& b);

kmer128_t operator& (kmer128_t const& a, kmer128_t const& b);

kmer128_t operator| (kmer128_t const& a, kmer128_t const& b);

kmer128_t operator^ (kmer128_t const& a, kmer128_t const& b);

kmer128_t operator- (kmer128_t const& a, int b);

kmer128_t operator<< (kmer128_t const& val, unsigned int shift);

kmer128_t operator>> (kmer128_t const& val, unsigned int shift);

//------------------------------------------------------------------------------------------------------------------

namespace debug {

struct triplet_t {
    uint64_t first, second, third;
};

// static uint64_t char_to_uint64(char c) { return (c >> 1) & 3; }

template <typename KMerType>
static KMerType char_to_uint(char c) {
    return constants::seq_nt4_table[static_cast<uint8_t>(c)] & 3;
}

[[maybe_unused]] static uint64_t extract_minimizer(uint64_t v) {return v;}
[[maybe_unused]] static uint64_t extract_minimizer(kmer128_t v) {return v.lower;}

template <typename KMerType>
[[maybe_unused]] static KMerType string_to_integer_no_reverse(const char* str, uint64_t k) 
{
    KMerType y;
    assert(k <= 64);
    y = 0;
    for (uint64_t i = 0; i != k; ++i) {
        // x += char_to_uint64(str[i]) << (2 * i);
        y = (y << 2) | char_to_uint<KMerType>(str[i]);
        // std::cerr << x << " " << y << "\n";
        // assert(x == y);
    }
    return y;
}

template <typename Hasher = hash64, typename KMerType>
static void print_hashes(std::string contig, uint64_t m, uint64_t seed) 
{
    for (uint64_t i = 0; i < contig.length() - m + 1; ++i) {
        KMerType pmmer = string_to_integer_no_reverse<KMerType>(&contig.data()[i], m);
        uint64_t hash = Hasher::hash(pmmer, seed).first();
        std::cerr << "[" << i << "] : " << pmmer << " " << hash << "\n";
    }
}

template <typename Hasher = hash64, typename KMerType>
static triplet_t compute_minimizer_triplet(KMerType kmer, uint64_t k, uint64_t m, uint64_t seed) 
{
    assert(m <= 32);
    assert(m <= k);
    uint64_t min_hash = uint64_t(-1);
    uint64_t minimizer = uint64_t(-1);
    KMerType mask = (static_cast<KMerType>(1) << (2 * m)) - 1;
    uint64_t pos = 0;
    for (uint64_t i = 0; i != k - m + 1; ++i) {
        uint64_t mmer = extract_minimizer(kmer & mask);
        uint64_t hash = Hasher::hash(mmer, seed).first();
        if (hash <= min_hash) {  // <= because during construction we take the left-most minimum. 
                                 //Here we start looking from the right so we need to update everytime.
            min_hash = hash;
            minimizer = mmer;
            pos = i;
        }
        kmer = kmer >> 2;
    }
    return triplet_t{minimizer, min_hash, k - (pos + m)};
}

template <typename Hasher = hash64, typename KMerType>
static void compute_minimizers_naive(std::string const& contig, uint64_t k, uint64_t m, uint64_t seed) 
{
    triplet_t prev = {uint64_t(-1), uint64_t(-1), uint64_t(-1)};
    for (std::size_t i = 0; i < contig.size() - k + 1; ++i) {
        KMerType integer_kmer = string_to_integer_no_reverse<KMerType>(&contig.data()[i], k);
        auto curr = compute_minimizer_triplet(integer_kmer, k, m, seed);
        if (prev.first != curr.first or (prev.third - curr.third > 1))
            std::cerr << curr.first << " " << curr.second << " " << curr.third << "\n";
        prev = curr;
    }
}

}  // namespace debug

}  // namespace lphash

#endif