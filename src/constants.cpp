#include "../include/constants.hpp"

namespace lphash {

namespace constants {

const std::array<uint8_t, 256> seq_nt4_table = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

}

bool operator< (mm_triplet_t const& a, mm_triplet_t const& b)
{
    return a.itself < b.itself;
}

std::ostream &operator<<(std::ostream &os, kmer128_t const& val) 
{
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

}