#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <array>
#include <string>
#include <iostream>

#include "../external/pthash/include/pthash.hpp"

namespace lphash {

struct kmer128_t {
    kmer128_t() {upper = 0; lower = 0;};
    kmer128_t(uint8_t v) {upper = 0; lower = static_cast<uint64_t>(v);};
    kmer128_t(int v) {upper = 0; lower = static_cast<uint64_t>(v);};
    kmer128_t(uint64_t v) {upper = 0; lower = v;};
    uint64_t upper;
    uint64_t lower;
};

#include "compile_constants.tpd"

namespace constants {

template <typename T> struct MaxKChooser;

template <> struct MaxKChooser<uint64_t>
{
    static uint64_t const value = 31;
};

template <> struct MaxKChooser<kmer128_t>
{
    static uint64_t const value = 63;
};

constexpr uint64_t max_k = MaxKChooser<kmer_t>::value;  // max *odd* size that can be packed into the given k-mer type
constexpr uint64_t seed = 1;
constexpr double c = 3.0;  // for PTHash
static const std::string default_tmp_dirname(".");
extern const std::array<uint8_t, 256> seq_nt4_table;

} // namespace constants

#pragma pack(push, 2)
struct mm_record_t {
    uint64_t itself;
    uint64_t id;
    uint8_t p1;
    uint8_t size;
};
#pragma pack(pop)

#pragma pack(push, 2)
struct mm_triplet_t {
    uint64_t itself;
    uint8_t p1;
    uint8_t size;
};
#pragma pack(pop)

struct hash64 : public pthash::murmurhash2_64 {
    static inline pthash::hash64 hash(kmer128_t val, uint64_t seed) {
        return pthash::MurmurHash2_64(reinterpret_cast<char const*>(&val), sizeof(val), seed);
    }
    static inline pthash::hash64 hash(pthash::byte_range range, uint64_t seed) {
        return murmurhash2_64::hash(range, seed);
    }
    static inline pthash::hash64 hash(std::string const& val, uint64_t seed) {
        return murmurhash2_64::hash(val, seed);
    }
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

std::string get_group_id();

template <class T>
void explicit_garbage_collect([[maybe_unused]] T obj) {}

//------------------------------------------------------------------------------------------------------------------

namespace debug {

struct triplet_t {
    uint64_t first, second, third;
};

static kmer_t char_to_uint(char c) {
    return constants::seq_nt4_table[static_cast<uint8_t>(c)] & 3;
}

[[maybe_unused]] static uint64_t extract_minimizer(uint64_t v) {return v;}
[[maybe_unused]] static uint64_t extract_minimizer(kmer128_t v) {return v.lower;}

[[maybe_unused]] static kmer_t string_to_integer_no_reverse(const char* str, uint64_t k) 
{
    kmer_t y;
    assert(k <= 64);
    y = 0;
    for (uint64_t i = 0; i != k; ++i) {
        y = (y << 2) | char_to_uint(str[i]);
    }
    return y;
}

template <typename Hasher = hash64>
static void print_hashes(std::string contig, uint64_t m, uint64_t seed) 
{
    for (uint64_t i = 0; i < contig.length() - m + 1; ++i) {
        kmer_t pmmer = string_to_integer_no_reverse(&contig.data()[i], m);
        uint64_t hash = Hasher::hash(pmmer, seed).first();
        std::cerr << "[" << i << "] : " << pmmer << " " << hash << "\n";
    }
}

template <typename Hasher = hash64>
static triplet_t compute_minimizer_triplet(kmer_t kmer, uint64_t k, uint64_t m, uint64_t seed) 
{
    assert(m <= 32);
    assert(m <= k);
    uint64_t min_hash = uint64_t(-1);
    uint64_t minimizer = uint64_t(-1);
    kmer_t mask = (static_cast<kmer_t>(1) << (2 * m)) - 1;
    uint64_t pos = 0;
    for (uint64_t i = 0; i != k - m + 1; ++i) {
        uint64_t mmer = extract_minimizer(kmer & mask);
        uint64_t hash = Hasher::hash(mmer, seed).first();
        if (hash <= min_hash) {  // <= because during construction we take the left-most minimum. 
                                 // Here we start looking from the right so we need to update every time.
            min_hash = hash;
            minimizer = mmer;
            pos = i;
        }
        kmer = kmer >> 2;
    }
    return triplet_t{minimizer, min_hash, k - (pos + m)};
}

template <typename Hasher = hash64>
static void compute_minimizers_naive(std::string const& contig, uint64_t k, uint64_t m, uint64_t seed) 
{
    triplet_t prev = {uint64_t(-1), uint64_t(-1), uint64_t(-1)};
    for (std::size_t i = 0; i < contig.size() - k + 1; ++i) {
        kmer_t integer_kmer = string_to_integer_no_reverse(&contig.data()[i], k);
        auto curr = compute_minimizer_triplet(integer_kmer, k, m, seed);
        if (prev.first != curr.first or (prev.third - curr.third > 1))
            std::cerr << curr.first << " " << curr.second << " " << curr.third << "\n";
        prev = curr;
    }
}

}  // namespace debug

}  // namespace lphash

#endif