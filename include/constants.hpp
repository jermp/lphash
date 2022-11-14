#pragma once

#include <array>
#include <string>
#include <iostream>

#include "../external/pthash/include/pthash.hpp"

namespace lphash {

struct kmer128_t {
    kmer128_t() : upper(0), lower(0) {}
    // kmer128_t(uint8_t v) : upper(0), lower(static_cast<uint64_t>(v)) {}
    // kmer128_t(int v) : upper(0), lower(static_cast<uint64_t>(v)) {}
    kmer128_t(uint64_t v) : upper(0), lower(v) {}
    uint64_t upper;
    uint64_t lower;
};

#include "compile_constants.tpd"

namespace constants {

template <typename T>
struct MaxKChooser;

template <>
struct MaxKChooser<uint64_t> {
    static uint64_t const value = 31;
};

template <>
struct MaxKChooser<kmer128_t> {
    static uint64_t const value = 63;
};

constexpr uint64_t max_k =
    MaxKChooser<kmer_t>::value;  // max *odd* size that can be packed into the given k-mer type
constexpr uint64_t default_pthash_seed = 1;
constexpr uint64_t default_seed = 42;
constexpr uint64_t default_num_threads = 1;
constexpr double c = 3.0;  // for PTHash
static const std::string default_tmp_dirname(".");
extern const std::array<uint8_t, 256> seq_nt4_table;

}  // namespace constants

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

bool operator<(mm_record_t const& a, mm_record_t const& b);
bool operator>(mm_record_t const& a, mm_record_t const& b);
bool operator<(mm_triplet_t const& a, mm_triplet_t const& b);
bool operator>(mm_triplet_t const& a, mm_triplet_t const& b);
bool operator==(kmer128_t const& a, kmer128_t const& b);
bool operator!=(kmer128_t const& a, kmer128_t const& b);
bool operator<(kmer128_t const& a, kmer128_t const& b);
bool operator>(kmer128_t const& a, kmer128_t const& b);
std::ostream& operator<<(std::ostream& os, kmer128_t const& val);
kmer128_t operator&(kmer128_t const& a, kmer128_t const& b);
kmer128_t operator|(kmer128_t const& a, kmer128_t const& b);
kmer128_t operator^(kmer128_t const& a, kmer128_t const& b);
kmer128_t operator-(kmer128_t const& a, int b);
kmer128_t operator<<(kmer128_t const& val, unsigned int shift);
kmer128_t operator>>(kmer128_t const& val, unsigned int shift);

//------------------------------------------------------------------------------------------------------------------

std::string get_group_id();

template <class T>
void explicit_garbage_collect([[maybe_unused]] T obj) {}

}  // namespace lphash
