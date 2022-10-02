#ifndef MM_CONTEXT_HPP
#define MM_CONTEXT_HPP

#include "constants.hpp"

namespace lphash {

namespace minimizer {

struct mm_quartet_t {
    mm_quartet_t() : hash(0), p1(0), size(0){};
    void clear() noexcept {
        itself = 0;  // AAA...A
        hash = std::numeric_limits<decltype(hash)>::max();
        p1 = std::numeric_limits<decltype(p1)>::max();
        size = std::numeric_limits<decltype(size)>::max();
    };

    friend std::ostream& operator<<(std::ostream& os, mm_quartet_t const& other) {
        return os << other.itself << " " << other.hash << ' ' << other.p1 << ' ' << other.size;
    }

    uint64_t hash;    // minimizer hash
    uint64_t itself;  // 2-bit minimizer itself
    uint32_t p1;      // position inside first k-mer of the super-k-mer
    uint32_t size;    // size (number of k-mers) in the super-k-mer
};

} // namespace minimizer

} // namespace lphash

#endif // MM_CONTEXT_HPP