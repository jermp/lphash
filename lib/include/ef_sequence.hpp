#pragma once

#include "../../external/pthash/include/encoders/bit_vector.hpp"
#include "../../external/pthash/include/encoders/compact_vector.hpp"
#include "../../external/pthash/include/encoders/darray.hpp"

namespace lphash {

template <typename Iterator>
class cumulative_iterator {
public:
    typedef std::size_t value_type;
    cumulative_iterator(Iterator& other) : inc(false), itr(other), cumulative_sum(*other) {}
    void operator++() {
        if (inc) cumulative_sum += *itr;
        ++itr;
        inc = true;
    }
    std::size_t operator*() {
        if (inc) {
            cumulative_sum += *itr;
            inc = false;
        }
        return cumulative_sum;
    }

private:
    bool inc;
    Iterator& itr;
    std::size_t cumulative_sum;
};

struct ef_sequence {
    ef_sequence() {}

    template <typename Iterator>
    void encode(Iterator begin, uint64_t n, uint64_t u) {
        if (n == 0) return;

        // if constexpr (encode_prefix_sum) {
        n = n + 1;  // because I will add a zero at the beginning
        // }

        uint64_t l = uint64_t((n && u / n) ? pthash::util::msb(u / n) : 0);
        pthash::bit_vector_builder bvb_high_bits(n + (u >> l) + 1);
        pthash::compact_vector::builder cv_builder_low_bits(n, l);

        uint64_t low_mask = (uint64_t(1) << l) - 1;
        uint64_t last = 0;
        // I add a zero at the beginning
        // if constexpr (encode_prefix_sum) {
        if (l) cv_builder_low_bits.push_back(0);
        bvb_high_bits.set(0, 1);
        n = n - 1;  // restore n
        // }
        for (size_t i = 0; i < n; ++i) {
            auto v = *begin;
            // if constexpr (encode_prefix_sum) {
            //     v = v + last;             // prefix sum
            if (i and v < last) {  // check the order
                std::cerr << "error at " << i << "/" << n << ":\n";
                std::cerr << "last " << last << "\n";
                std::cerr << "current " << static_cast<uint64_t>(v) << "\n";
                throw std::runtime_error("ef_sequence is not sorted");
            }
            if (l) cv_builder_low_bits.push_back(v & low_mask);
            bvb_high_bits.set((v >> l) + i + true, 1);  // encode_prefix_sum, 1);
            last = v;
            ++begin;
        }

        pthash::bit_vector(&bvb_high_bits).swap(m_high_bits);
        cv_builder_low_bits.build(m_low_bits);
        m_high_bits_d1.build(m_high_bits);
    }

    inline uint64_t access(uint64_t i) const {
        assert(i < size());
        return ((m_high_bits_d1.select(m_high_bits, i) - i) << m_low_bits.width()) |
               m_low_bits.access(i);
    }

    inline std::pair<uint64_t, uint64_t> pair(uint64_t i) const {
        assert(i < size());  // and encode_prefix_sum);
        uint64_t low1 = m_low_bits.access(i);
        uint64_t low2 = m_low_bits.access(i + 1);
        uint64_t l = m_low_bits.width();
        uint64_t pos = m_high_bits_d1.select(m_high_bits, i);
        uint64_t h1 = pos - i;
        uint64_t h2 = pthash::bit_vector::unary_iterator(m_high_bits, pos + 1).next() - i - 1;
        uint64_t val1 = (h1 << l) | low1;
        uint64_t val2 = (h2 << l) | low2;
        return {val1, val2};
    }

    inline uint64_t diff(uint64_t i) const {
        auto [val1, val2] = pair(i);
        return val2 - val1;
    }

    inline uint64_t size() const { return m_low_bits.size(); }

    uint64_t num_bits() const {
        return 8 * (m_high_bits.bytes() + m_high_bits_d1.bytes() + m_low_bits.bytes());
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_high_bits);
        visitor.visit(m_high_bits_d1);
        visitor.visit(m_low_bits);
    }

private:
    pthash::bit_vector m_high_bits;
    pthash::darray1 m_high_bits_d1;
    pthash::compact_vector m_low_bits;
};

}  // namespace lphash