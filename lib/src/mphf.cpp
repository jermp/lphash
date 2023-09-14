#include "../include/mphf.hpp"

namespace lphash {
namespace mphf {

interface::interface() : k(0), m(0), mm_seed(0), nkmers(0), distinct_minimizers(0), max_ram(0) {}

uint8_t interface::get_k() const noexcept {
    return k;
}

uint8_t interface::get_m() const noexcept {
    return m;
}

uint64_t interface::get_minimizer_L0() const noexcept { return distinct_minimizers; }

uint64_t interface::get_kmer_count() const noexcept { return nkmers; }

uint64_t interface::get_minimizer_order(uint64_t mm) const { return minimizer_order(mm); }

void interface::build_minimizers_mphf(
    external_memory_vector<mm_triplet_t, false>::const_iterator& mm_itr,
    std::size_t number_of_distinct_minimizers) {
    mm_itr_t dummy_itr(mm_itr);
    distinct_minimizers = number_of_distinct_minimizers;
    minimizer_order.build_in_external_memory(std::move(dummy_itr), distinct_minimizers, mphf_configuration);
}

void interface::build_fallback_mphf(external_memory_vector<kmer_t, false>::const_iterator& km_itr, std::size_t number_of_colliding_kmers) {
    km_itr_t dummy_itr(km_itr);
    fallback_kmer_order.build_in_external_memory(dummy_itr, number_of_colliding_kmers, mphf_configuration);
}

mm_itr_t::mm_itr_t(external_memory_vector<mm_triplet_t, false>::const_iterator& mm_itr)
    : m_iterator(mm_itr) {}

uint64_t mm_itr_t::operator*() const { return (*m_iterator).itself; }

void mm_itr_t::operator++() { ++m_iterator; }

km_itr_t::km_itr_t(external_memory_vector<kmer_t, false>::const_iterator& km_itr)
    : m_iterator(km_itr) {}

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

} // namespace mphf
}  // namespace lphash