#include "../include/mphf_utils.hpp"

namespace lphash {

mm_itr_t::mm_itr_t(external_memory_vector<mm_triplet_t>::const_iterator& mm_itr)
    : m_iterator(mm_itr) {}

uint64_t mm_itr_t::operator*() const { return (*m_iterator).itself; }

void mm_itr_t::operator++() { ++m_iterator; }

km_itr_t::km_itr_t(external_memory_vector<kmer_t>::const_iterator& km_itr) : m_iterator(km_itr) {}

kmer_t const& km_itr_t::operator*() const { return (*m_iterator); }

void km_itr_t::operator++() { ++m_iterator; }

pos_itr_t::pos_itr_t(external_memory_vector<mm_triplet_t>::const_iterator& km_itr)
    : m_iterator(km_itr) {}

uint8_t const& pos_itr_t::operator*() const { return (*m_iterator).p1; }

void pos_itr_t::operator++() { ++m_iterator; }

size_itr_t::size_itr_t(external_memory_vector<mm_triplet_t>::const_iterator& km_itr)
    : m_iterator(km_itr) {}

uint8_t const& size_itr_t::operator*() const { return (*m_iterator).size; }

void size_itr_t::operator++() { ++m_iterator; }

}  // namespace lphash