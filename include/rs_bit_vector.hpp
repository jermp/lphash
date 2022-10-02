#pragma once

#include <algorithm>
#include <vector>

#include "../external/pthash/external/essentials/include/essentials.hpp"
#include "../external/pthash/include/encoders/bit_vector.hpp"

class rs_bit_vector : public pthash::bit_vector {
public:
    rs_bit_vector() : pthash::bit_vector() {}

    void build(pthash::bit_vector_builder* bvb, bool with_select_hints) {
        pthash::bit_vector::build(bvb);
        build_indices(with_select_hints);
    }

    rs_bit_vector(pthash::bit_vector_builder* bvb, bool with_select_hints) {
        build(bvb, with_select_hints);
    }

    // void swap(rs_bit_vector& other) {
    //     pthash::bit_vector::swap(other);
    //     m_block_rank_pairs.swap(other.m_block_rank_pairs);
    //     m_select_hints.swap(other.m_select_hints);
    // }

    inline uint64_t num_ones() const { return *(m_block_rank_pairs.end() - 2); }

    inline uint64_t num_zeros() const { return size() - num_ones(); }

    inline uint64_t rank(uint64_t pos) const {
        assert(pos <= size());
        if (pos == size()) { return num_ones(); }

        uint64_t sub_block = pos / 64;
        uint64_t r = sub_block_rank(sub_block);
        uint64_t sub_left = pos % 64;
        if (sub_left) { r += pthash::util::popcount(m_bits[sub_block] << (64 - sub_left)); }
        return r;
    }

    inline uint64_t rank0(uint64_t pos) const { return pos - rank(pos); }

    inline uint64_t select(uint64_t n) const {
        assert(n < num_ones());
        uint64_t a = 0;
        uint64_t b = num_blocks();
        if (m_select_hints.size()) {
            uint64_t chunk = n / select_ones_per_hint;
            if (chunk != 0) a = m_select_hints[chunk - 1];
            b = m_select_hints[chunk] + 1;
        }

        uint64_t block = 0;
        while (b - a > 1) {
            uint64_t mid = a + (b - a) / 2;
            uint64_t x = block_rank(mid);
            if (x <= n) {
                a = mid;
            } else {
                b = mid;
            }
        }
        block = a;

        assert(block < num_blocks());
        uint64_t block_offset = block * block_size;
        uint64_t cur_rank = block_rank(block);
        assert(cur_rank <= n);

        uint64_t rank_in_block_parallel = (n - cur_rank) * ones_step_9;
        uint64_t sub_ranks = sub_block_ranks(block);
        uint64_t sub_block_offset =
            uleq_step_9(sub_ranks, rank_in_block_parallel) * ones_step_9 >> 54 & 0x7;
        cur_rank += sub_ranks >> (7 - sub_block_offset) * 9 & 0x1FF;
        assert(cur_rank <= n);

        uint64_t word_offset = block_offset + sub_block_offset;
        return word_offset * 64 + pthash::util::select_in_word(m_bits[word_offset], n - cur_rank);
    }

    static const uint64_t ones_step_9 =
        1ULL << 0 | 1ULL << 9 | 1ULL << 18 | 1ULL << 27 | 1ULL << 36 | 1ULL << 45 | 1ULL << 54;
    static const uint64_t msbs_step_9 = 0x100ULL * ones_step_9;
    inline static uint64_t uleq_step_9(uint64_t x, uint64_t y) {
        return (((((y | msbs_step_9) - (x & ~msbs_step_9)) | (x ^ y)) ^ (x & ~y)) & msbs_step_9) >>
               8;
    }

    uint64_t bytes() const {
        return pthash::bit_vector::bytes() + essentials::vec_bytes(m_block_rank_pairs) +
               essentials::vec_bytes(m_select_hints);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        pthash::bit_vector::visit(visitor);
        visitor.visit(m_block_rank_pairs);
        visitor.visit(m_select_hints);
    }

protected:
    inline uint64_t num_blocks() const { return m_block_rank_pairs.size() / 2 - 1; }

    inline uint64_t block_rank(uint64_t block) const { return m_block_rank_pairs[block * 2]; }

    inline uint64_t sub_block_rank(uint64_t sub_block) const {
        uint64_t r = 0;
        uint64_t block = sub_block / block_size;
        r += block_rank(block);
        uint64_t left = sub_block % block_size;
        r += sub_block_ranks(block) >> ((7 - left) * 9) & 0x1FF;
        return r;
    }

    inline uint64_t sub_block_ranks(uint64_t block) const {
        return m_block_rank_pairs[block * 2 + 1];
    }

    inline uint64_t block_rank0(uint64_t block) const {
        return block * block_size * 64 - m_block_rank_pairs[block * 2];
    }

    void build_indices(bool with_select_hints) {
        {
            std::vector<uint64_t> block_rank_pairs;
            uint64_t next_rank = 0;
            uint64_t cur_subrank = 0;
            uint64_t subranks = 0;
            block_rank_pairs.push_back(0);
            for (uint64_t i = 0; i < m_bits.size(); ++i) {
                uint64_t word_pop = pthash::util::popcount(m_bits[i]);
                uint64_t shift = i % block_size;
                if (shift) {
                    subranks <<= 9;
                    subranks |= cur_subrank;
                }
                next_rank += word_pop;
                cur_subrank += word_pop;

                if (shift == block_size - 1) {
                    block_rank_pairs.push_back(subranks);
                    block_rank_pairs.push_back(next_rank);
                    subranks = 0;
                    cur_subrank = 0;
                }
            }
            uint64_t left = block_size - m_bits.size() % block_size;
            for (uint64_t i = 0; i < left; ++i) {
                subranks <<= 9;
                subranks |= cur_subrank;
            }
            block_rank_pairs.push_back(subranks);

            if (m_bits.size() % block_size) {
                block_rank_pairs.push_back(next_rank);
                block_rank_pairs.push_back(0);
            }

            m_block_rank_pairs.swap(block_rank_pairs);
        }

        if (with_select_hints) {
            std::vector<uint64_t> select_hints;
            uint64_t cur_ones_threshold = select_ones_per_hint;
            for (uint64_t i = 0; i < num_blocks(); ++i) {
                if (block_rank(i + 1) > cur_ones_threshold) {
                    select_hints.push_back(i);
                    cur_ones_threshold += select_ones_per_hint;
                }
            }
            select_hints.push_back(num_blocks());
            m_select_hints.swap(select_hints);
        }
    }

    static const uint64_t block_size = 8;                              // in 64bit words
    static const uint64_t select_ones_per_hint = 64 * block_size * 2;  // must be > block_size * 64
    static const uint64_t select_zeros_per_hint = select_ones_per_hint;

    std::vector<uint64_t> m_block_rank_pairs;
    std::vector<uint64_t> m_select_hints;
};