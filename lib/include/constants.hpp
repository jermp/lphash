#pragma once

#include <array>
#include <iostream>

#include "../../external/pthash/include/pthash.hpp"

namespace lphash {

#include "../../main/compile_constants.tpd" // how to fix this?

namespace constants {

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

// struct fallback_hasher {
//     typedef pthash::hash128 hash_type;
//     static inline pthash::hash128 hash(kmer_t val, uint64_t seed) {
//         return {pthash::MurmurHash2_64(reinterpret_cast<char const*>(&val), sizeof(val), seed),
//                 pthash::MurmurHash2_64(reinterpret_cast<char const*>(&val), sizeof(val), ~seed)};
//     }
// };

struct fallback_hasher {
    typedef pthash::hash64 hash_type;
    static inline pthash::hash64 hash(kmer_t val, uint64_t seed) {
        if constexpr (sizeof(val) == sizeof(uint64_t)) {
            return pthash::MurmurHash2_64(reinterpret_cast<char const*>(&val), sizeof(val), seed);
        } else {
            uint64_t low = static_cast<uint64_t>(val);
            uint64_t high = static_cast<uint64_t>(val >> 64);
            return pthash::MurmurHash2_64(reinterpret_cast<char const*>(&low), sizeof(uint64_t),
                                          seed) ^
                   pthash::MurmurHash2_64(reinterpret_cast<char const*>(&high), sizeof(uint64_t),
                                          ~seed);
        }
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
