#pragma once

#include <array>
#include <string>
#include <iostream>

#include "../external/pthash/include/pthash.hpp"

namespace lphash {

#include "compile_constants.tpd"

namespace constants {

/* max *odd* size that can be packed into the given k-mer type */
constexpr uint64_t max_k = (sizeof(kmer_t) * 8) / 2 - 1;
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

static bool operator<(mm_record_t const& a, mm_record_t const& b) { return a.itself < b.itself; }
static bool operator>(mm_record_t const& a, mm_record_t const& b) { return a.itself > b.itself; }
static bool operator<(mm_triplet_t const& a, mm_triplet_t const& b) { return a.itself < b.itself; }
static bool operator>(mm_triplet_t const& a, mm_triplet_t const& b) { return a.itself < b.itself; }

struct fallback_hasher {
    typedef pthash::hash128 hash_type;
    static inline pthash::hash128 hash(kmer_t val, uint64_t seed) {
        return {pthash::MurmurHash2_64(reinterpret_cast<char const*>(&val), sizeof(val), seed),
                pthash::MurmurHash2_64(reinterpret_cast<char const*>(&val), sizeof(val), ~seed)};
    }
};

typedef pthash::single_phf<pthash::murmurhash2_64, pthash::dictionary_dictionary, true>
    pthash_minimizers_mphf_t;
typedef pthash::single_phf<fallback_hasher, pthash::dictionary_dictionary, true>
    pthash_fallback_mphf_t;

static std::string get_group_id() {
    return std::to_string(pthash::clock_type::now().time_since_epoch().count());
}

template <class T>
void explicit_garbage_collect([[maybe_unused]] T obj) {}

}  // namespace lphash
