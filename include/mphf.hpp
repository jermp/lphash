#ifndef MPHF_HPP
#define MPHF_HPP

#include "../include/constants.hpp"
#include "quartet_wtree.hpp"

namespace lphash {

class mphf {
    public:
        mphf();
        mphf(uint8_t klen, uint8_t mm_size, uint64_t seed, uint64_t total_number_of_kmers, uint8_t nthreads, std::string temporary_directory = "", bool verbose = false);
        void build_minimizers_mphf(std::vector<mm_triplet_t>& minimizers);
        template<typename KMerType> void build_fallback_mphf(std::vector<KMerType>& colliding_kmers);
        std::vector<uint64_t> build_inverted_index(std::vector<mm_triplet_t>& minimizers);
        uint64_t get_minimizer_L0() const noexcept;
        uint64_t get_kmer_count() const noexcept;
        std::vector<uint64_t> operator() (std::string const& contig) const;
        void print_statistics() const noexcept;
        template <typename Visitor> void visit(Visitor& visitor);
        pthash_mphf_t minimizer_order;
        pthash_mphf_t fallback_kmer_order;
        
    private:
        pthash::build_configuration mphf_configuration;
        uint8_t k, m;
        uint64_t mm_seed;
        uint64_t nkmers;
        uint64_t distinct_minimizers;
        uint64_t n_maximal, right_coll_sizes_start, none_sizes_start, none_pos_start; // Left positions | right + coll sizes | none sizes | none positions
        quartet_wtree wtree;
        pthash::ef_sequence<true> sizes_and_positions;
};

template <typename KMerType>
void mphf::build_fallback_mphf(std::vector<KMerType>& colliding_kmers)
{
    std::sort(colliding_kmers.begin(), colliding_kmers.end());
    auto it = std::unique(colliding_kmers.begin(), colliding_kmers.end());
    {
        std::ofstream kmers("kmers.txt");
        for (auto kmer : colliding_kmers) {
            kmers << kmer << "\n";
        }
    }
    assert(it == colliding_kmers.end()); // if false: there are some duplicates in unbucketable_kmers
    fallback_kmer_order.build_in_external_memory(colliding_kmers.begin(), colliding_kmers.size(), mphf_configuration);
}

template <typename Visitor>
void mphf::visit(Visitor& visitor) {
    //visitor.visit(mphf_configuration);
    visitor.visit(k);
    visitor.visit(m);
    visitor.visit(mm_seed);
    visitor.visit(nkmers);
    visitor.visit(distinct_minimizers);
    visitor.visit(n_maximal);
    visitor.visit(right_coll_sizes_start);
    visitor.visit(none_sizes_start);
    visitor.visit(none_pos_start);
    visitor.visit(wtree);
    visitor.visit(sizes_and_positions);
}

bool check_collisions(mphf const& hf, std::string const& contig, pthash::bit_vector_builder& population);
bool check_perfection(mphf const& hf, pthash::bit_vector_builder& population);

} // namespace lphash

#endif // MPHF_HPP