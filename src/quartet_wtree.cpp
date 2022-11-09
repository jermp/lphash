#include "../include/quartet_wtree.hpp"

namespace lphash {

quartet_wtree_builder::quartet_wtree_builder(std::size_t size)
    : m_root_size(0), m_left_right_size(0), m_max_none_size(0) {
    root.reserve(size);
    left_right.reserve(size);
    max_none.reserve(size);
}

void quartet_wtree_builder::push_back(MinimizerType type) {
    switch (type) {
        case MAXIMAL:
            root.push_back(true);
            max_none.push_back(false);
            ++m_max_none_size;
            ++m_root_size;
            break;
        case LEFT:
            root.push_back(false);
            left_right.push_back(false);
            ++m_left_right_size;
            ++m_root_size;
            break;
        case RIGHT_OR_COLLISION:
            root.push_back(false);
            left_right.push_back(true);
            ++m_left_right_size;
            ++m_root_size;
            break;
        case NONE:
            root.push_back(true);
            max_none.push_back(true);
            ++m_max_none_size;
            ++m_root_size;
            break;
        default:
            break;
    }
}

void quartet_wtree_builder::freeze() {
    root.resize(m_root_size);
    left_right.resize(m_left_right_size);
    max_none.resize(m_max_none_size);
}

void quartet_wtree::build(quartet_wtree_builder& unfrozen) {
    unfrozen.freeze();
    root.build(&unfrozen.root, false);  // false?
    left_right.build(&unfrozen.left_right, false);
    max_none.build(&unfrozen.max_none, false);
}

MinimizerType quartet_wtree::operator[](uint64_t idx) const {
    bool msb = root[idx];
    auto r = rank_switch(msb, root, idx);
    bool lsb;
    if (msb)
        lsb = max_none[r];
    else
        lsb = left_right[r];
    return static_cast<MinimizerType>((static_cast<uint64_t>(msb) << 1) |
                                      static_cast<uint64_t>(lsb));
}

std::size_t quartet_wtree::rank(MinimizerType type, std::size_t idx) const {
    bool sw_leaf = static_cast<uint64_t>(type) & 1;
    bool sw_root = static_cast<uint64_t>(type) >> 1;
    auto leaf_idx = rank_switch(sw_root, root, idx);
    if (leaf_idx == 0) return 0;
    --leaf_idx;
    if (sw_root)
        leaf_idx = rank_switch(sw_leaf, max_none, leaf_idx);
    else
        leaf_idx = rank_switch(sw_leaf, left_right, leaf_idx);
    if (leaf_idx == 0)
        return 0;
    else
        return leaf_idx;
}

std::pair<MinimizerType, std::size_t> quartet_wtree::rank_of(std::size_t idx) const {
    bool msb = root[idx];
    bool lsb;
    std::pair<MinimizerType, std::size_t> toRet;
    auto r = rank_switch(msb, root, idx);
    if (msb) {
        lsb = max_none[r];
        toRet.second = rank_switch(lsb, max_none, r);
    } else {
        lsb = left_right[r];
        toRet.second = rank_switch(lsb, left_right, r);
    }
    toRet.first =
        static_cast<MinimizerType>((static_cast<uint64_t>(msb) << 1) | static_cast<uint64_t>(lsb));
    return toRet;
}

std::size_t quartet_wtree::rank_switch(bool type, rs_bit_vector const& vec, std::size_t idx) const {
    if (type)
        return vec.rank(idx);
    else
        return vec.rank0(idx);
}

std::size_t quartet_wtree::num_bits() const {
    return (root.bytes() + left_right.bytes() + max_none.bytes()) * 8;
}

}  // namespace lphash