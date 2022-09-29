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
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

}

struct kmer128_t {
    kmer128_t() {upper = 0; lower = 0;};
    kmer128_t(uint8_t v) {upper = 0; lower = static_cast<uint64_t>(v);};
    kmer128_t(int v) {upper = 0; lower = static_cast<uint64_t>(v);};
    kmer128_t(uint64_t v) {upper = 0; lower = v;};
    // operator uint64_t() const {return lower;};
    uint64_t upper;
    uint64_t lower;
};

std::ostream &operator<<(std::ostream &os, kmer128_t const& val) {
    return os << "(" << val.upper << ", " << val.lower << ")";
}

bool operator< (kmer128_t const& a, kmer128_t const& b)
{
    if (a.upper == b.upper) return a.lower < b.lower;
    return (a.upper < b.upper);
}

bool operator== (kmer128_t const& a, kmer128_t const& b)
{
    return a.upper == b.upper && a.lower == b.lower;
}

kmer128_t operator& (kmer128_t const& a, kmer128_t const& b)
{
    kmer128_t res;
    res.upper = a.upper & b.upper;
    res.lower = a.lower & b.lower;
    return res;
}

kmer128_t operator| (kmer128_t const& a, kmer128_t const& b)
{
    kmer128_t res;
    res.upper = a.upper | b.upper;
    res.lower = a.lower | b.lower;
    return res;
}

kmer128_t operator^ (kmer128_t const& a, kmer128_t const& b)
{
    kmer128_t res;
    res.upper = a.upper ^ b.upper;
    res.lower = a.lower ^ b.lower;
    return res;
}

kmer128_t operator- (kmer128_t const& a, int b)
{
    kmer128_t res;
    res.lower = a.lower - b;
    if (res.lower > a.lower) res.upper = a.upper - 1;
    return res;
}

kmer128_t operator<< (kmer128_t const& val, unsigned int shift) 
{
    kmer128_t res;;
    if (shift < 64) {
        uint64_t mask = ~((1ULL << shift) - 1);
        res.upper = (val .lower & mask) >> (64 - shift);
        res.lower = val.lower << shift;
        res.upper |= val.upper << shift;
    } else {
        res.lower = 0;
        res.upper = val.lower << (shift % 64);
    }
    return res;
}

kmer128_t operator>> (kmer128_t const& val, unsigned int shift) 
{
    kmer128_t res;
    if (shift < 64) {
        uint64_t mask = (1ULL << shift) - 1;
        res.lower = (val.upper & mask) << (64 - shift);
        res.lower |= val.lower >> shift;
        res.upper = val.upper >> shift;
    } else {
        res.lower = val.upper >> (shift % 64);
        res.upper = 0;
    }
    return res;
}

typedef pthash::single_phf<pthash::murmurhash2_64,         // base hasher
                           pthash::dictionary_dictionary,  // encoder type
                           true                            // minimal output
                           >
    pthash_mphf_type;

// struct murmurhash2_64 { //my hash function for hashing minimizers
//     // typedef pthash::hash64 hash_type;
//     static inline uint64_t /*pthash::hash64*/ hash(uint64_t val, uint64_t seed) {
//         return pthash::MurmurHash2_64(reinterpret_cast<char const*>(&val), sizeof(val), seed);
//     }
//     // static inline uint64_t /*pthash::hash64*/ hash(kmer128_t val, uint64_t seed) {
//     //     return pthash::MurmurHash2_64(reinterpret_cast<char const*>(&val), sizeof(val), seed);
//     // }
// };

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

typedef pthash::single_phf<hash64, pthash::dictionary_dictionary, true> lphash_mphf_type;

typedef kmer128_t kmer_t;
// typedef uint64_t kmer_t;

// struct murmurhash2_128 {
//     typedef pthash::hash128 hash_type;
//     static inline pthash::hash128 hash(__uint128_t val, uint64_t seed) {
//         return {pthash::MurmurHash2_64(reinterpret_cast<char const*>(&val), 8, seed),
//                 pthash::MurmurHash2_64(reinterpret_cast<char const*>(&val) + 8, 8, ~seed)};
//     }
//     // static inline __uint128_t hash(__uint128_t val, uint64_t seed) {
//     //     uint8_t const* dummy = reinterpret_cast<uint8_t const*>(&val);
//     //     pthash::byte_range vrng;
//     //     vrng.begin = dummy;
//     //     vrng.end = &dummy[sizeof(__uint128_t) - 1];
//     //     pthash::murmurhash2_128 hasher;
//     //     pthash::hash128 hval = hasher.hash(vrng, seed);
//     //     __uint128_t to_ret = hval.first();
//     //     to_ret = (to_ret << 64) | hval.second();
//     //     return to_ret;
//     // }
// };

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

template <typename KMerType>
static KMerType char_to_uint(char c) {
    return constants::seq_nt4_table[static_cast<uint8_t>(c)] & 3;
}

[[maybe_unused]] static uint64_t extract_minimizer(uint64_t v) {return v;}
[[maybe_unused]] static uint64_t extract_minimizer(kmer128_t v) {return v.lower;}

template <typename KMerType>
[[maybe_unused]] static KMerType string_to_integer_no_reverse(const char const* str, uint64_t k) 
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
        // std::cerr << std::bitset<64>(uint64_kmer) << "\n";
        auto curr = compute_minimizer_triplet(integer_kmer, k, m, seed);
        if (prev.first != curr.first or (prev.third - curr.third > 1))
            std::cerr << curr.first << " " << curr.second << " " << curr.third << "\n";
        prev = curr;
    }
}

}  // namespace debug

}  // namespace lphash

#endif