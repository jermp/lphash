#ifndef QUARTET_WTREE
#define QUARTET_WTREE

#include "../include/rs_bit_vector.hpp"

namespace lphash {

enum MinimizerType {LEFT, RIGHT_OR_COLLISION, MAXIMAL, NONE};

class quartet_wtree;

class quartet_wtree_builder
{
    public:
        quartet_wtree_builder(std::size_t initial_size);
        void push_back(MinimizerType type);
        void freeze();

    private:
        pthash::bit_vector_builder root;
        pthash::bit_vector_builder left_right;
        pthash::bit_vector_builder max_none;
        std::size_t m_root_size;
        std::size_t m_left_right_size;
        std::size_t m_max_none_size;
        friend quartet_wtree;
};

class quartet_wtree 
{
    public:
        quartet_wtree() {};
        void build(quartet_wtree_builder& unfrozen);
        MinimizerType operator[](uint64_t idx) const;
        std::size_t rank(MinimizerType type, std::size_t idx) const;
        std::pair<MinimizerType, std::size_t> rank_of(std::size_t idx) const;
    private:
        rs_bit_vector root, left_right, max_none;
        std::size_t rank_switch(bool type, rs_bit_vector const& vec, std::size_t idx) const;
};

}

#endif //QUARTET_WTREE