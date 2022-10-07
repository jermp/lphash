#include "../include/mphf.hpp"

#include "../include/prettyprint.hpp"

namespace lphash {

class vector_mm_triplet_to_pthash_itr_adapter : std::forward_iterator_tag {
public:
    typedef mm_triplet_t value_type;
    vector_mm_triplet_to_pthash_itr_adapter(std::vector<mm_triplet_t>::iterator begin,
                                            std::vector<mm_triplet_t>::iterator end)
        : begin(begin), end(end), current(begin){};
    inline uint64_t operator*() const { return current->itself; };
    inline void operator++() {
        uint64_t prev_mm = current->itself;
        while (current != end && current->itself == prev_mm) { ++current; }
    };

private:
    std::vector<mm_triplet_t>::iterator begin, end, current;
};

mphf::mphf()
    : k(0)
    , m(0)
    , mm_seed(0)
    , nkmers(0)
    , distinct_minimizers(0)
    , n_maximal(0)
    , right_coll_sizes_start(0)
    , none_sizes_start(0)
    , none_pos_start(0) {
    mphf_configuration.minimal_output = true;
    mphf_configuration.seed = constants::seed;
    mphf_configuration.c = constants::c;
    mphf_configuration.alpha = 0.94;
    mphf_configuration.verbose_output = false;
    mphf_configuration.num_threads = 0;
    mphf_configuration.tmp_dir = "";
};

mphf::mphf(uint8_t klen, uint8_t mm_size, uint64_t seed, uint64_t total_number_of_kmers,
           uint8_t nthreads, std::string temporary_directory, bool verbose)
    : k(klen)
    , m(mm_size)
    , mm_seed(seed)
    , nkmers(total_number_of_kmers)
    , distinct_minimizers(0)
    , n_maximal(0)
    , right_coll_sizes_start(0)
    , none_sizes_start(0)
    , none_pos_start(0) {
    mphf_configuration.minimal_output = true;
    mphf_configuration.seed = constants::seed;
    mphf_configuration.c = constants::c;
    mphf_configuration.alpha = 0.94;
    mphf_configuration.verbose_output = verbose;
    mphf_configuration.num_threads = nthreads;
    if (temporary_directory != "") {
        mphf_configuration.tmp_dir = temporary_directory;
        essentials::create_directory(temporary_directory);
    }
}

void mphf::build_minimizers_mphf(std::vector<mm_triplet_t>& minimizers) {
    std::sort(minimizers.begin(), minimizers.end());
    distinct_minimizers = 0;
    for (auto it = minimizers.begin(), prev = minimizers.begin(); it != minimizers.end(); ++it) {
        // FIXME get number of distinct minimizer -> rework pthash interface for this
        if (prev->itself != it->itself) {
            ++distinct_minimizers;
            prev = it;
        }
    }
    if (minimizers.size()) ++distinct_minimizers;
    auto begin = vector_mm_triplet_to_pthash_itr_adapter(minimizers.begin(), minimizers.end());
    minimizer_order.build_in_external_memory(begin, distinct_minimizers, mphf_configuration);
}

void mphf::build_fallback_mphf(std::vector<kmer_t>& colliding_kmers) {
    std::sort(colliding_kmers.begin(), colliding_kmers.end());
    [[maybe_unused]] auto it = std::unique(colliding_kmers.begin(), colliding_kmers.end());
    assert(it ==
           colliding_kmers.end());  // if false: there are some duplicates in unbucketable_kmers
    fallback_kmer_order.build_in_external_memory(colliding_kmers.begin(), colliding_kmers.size(),
                                                 mphf_configuration);
}

std::vector<uint64_t> mphf::build_inverted_index(std::vector<mm_triplet_t>& minimizers) {
    auto mphf_compare = [this](mm_triplet_t const& a, mm_triplet_t const& b) {
        return minimizer_order(a.itself) < minimizer_order(b.itself);
    };
    std::sort(minimizers.begin(), minimizers.end(), mphf_compare);

    bool duplicate = false;  // fix corner case of last minimizer without a successor
    n_maximal = 0;
    quartet_wtree_builder wtb(minimizers.size());
    std::vector<uint64_t> colliding_minimizers;
    std::vector<uint64_t> right_or_collision_sizes;
    std::vector<uint64_t> left_positions;
    std::vector<uint64_t> none_positions, none_sizes;
    auto insert_unique_minimizer = [&](std::size_t i) {
        if (minimizers[i].p1 == k - m) {
            if (minimizers[i].size == k - m + 1) {
                wtb.push_back(MAXIMAL);
                ++n_maximal;
            } else {
                wtb.push_back(RIGHT_OR_COLLISION);
                right_or_collision_sizes.push_back(minimizers[i].size);
            }
        } else {
            if (minimizers[i].p1 == minimizers[i].size - 1) {
                wtb.push_back(LEFT);
                left_positions.push_back(
                    minimizers[i].p1 +
                    1);  // +1 because with p1 == 0 we have 1 k-mer in the prefix sum
            } else {
                wtb.push_back(NONE);
                none_positions.push_back(minimizers[i].p1);
                none_sizes.push_back(minimizers[i].size);
            }
        }
    };
    if (minimizers.size() == 0) return colliding_minimizers;
    for (std::size_t i = 0; i < minimizers.size() - 1; ++i) {
        if (minimizers[i].itself != minimizers[i + 1].itself) {  // unique minimizer?
            insert_unique_minimizer(i);
            duplicate = false;
        } else {  // collision
            wtb.push_back(RIGHT_OR_COLLISION);
            colliding_minimizers.push_back(minimizers[i].itself);
            right_or_collision_sizes.push_back(0);
            for (std::size_t j = i + 1;
                 j < minimizers.size() && minimizers[j].itself == minimizers[i].itself; ++j) {
                ++i;
            }
            duplicate = true;
        }
    }
    if (!duplicate) insert_unique_minimizer(minimizers.size() - 1);
    assert(none_positions.size() == none_sizes.size());
    if (mphf_configuration.verbose_output) {
        double maximal = static_cast<double>(n_maximal) / minimizers.size() * 100;
        double left = static_cast<double>(left_positions.size()) / minimizers.size() * 100;
        double right =
            static_cast<double>(right_or_collision_sizes.size() - colliding_minimizers.size()) /
            minimizers.size() * 100;
        double none = static_cast<double>(none_positions.size()) / minimizers.size() * 100;
        double ambiguous =
            static_cast<double>(colliding_minimizers.size()) / minimizers.size() * 100;
        std::cerr << "Percentage of maximal super-k-mers: " << maximal << "%\n";
        std::cerr << "Percentage of left-maximal super-k-mers: " << left << "%\n";
        std::cerr << "Percentage of right-maximal super-k-mers : " << right << "%\n";
        std::cerr << "Percentage of unclassified super-k-mers: " << none << "%\n";
        std::cerr << "Percentage of ambiguous minimizers: " << ambiguous << "%\n";
    }
    wtree.build(wtb);

    right_coll_sizes_start = left_positions.size();
    left_positions.insert(left_positions.end(), right_or_collision_sizes.begin(),
                          right_or_collision_sizes.end());
    right_or_collision_sizes.clear();

    none_sizes_start = left_positions.size();
    left_positions.insert(left_positions.end(), none_sizes.begin(), none_sizes.end());
    none_sizes.clear();

    none_pos_start = left_positions.size();
    left_positions.insert(left_positions.end(), none_positions.begin(), none_positions.end());
    none_positions.clear();

    sizes_and_positions.encode(left_positions.begin(), left_positions.size());

    return colliding_minimizers;
}

uint64_t mphf::get_minimizer_L0() const noexcept { return distinct_minimizers; }

uint64_t mphf::get_kmer_count() const noexcept { return nkmers; }

uint64_t mphf::num_bits() const noexcept {
    auto mm_mphf_size_bits = minimizer_order.num_bits();
    auto triplet_tree_size_bits = wtree.num_bits();
    auto elias_sequence_size_bits = (sizeof(n_maximal) + sizeof(right_coll_sizes_start) +
                                     sizeof(none_sizes_start) + sizeof(none_pos_start)) *
                                        8 +
                                    sizes_and_positions.num_bits();
    auto kmer_mphf_size_bits = fallback_kmer_order.num_bits();
    auto total_bit_size = mm_mphf_size_bits + triplet_tree_size_bits + elias_sequence_size_bits +
                          kmer_mphf_size_bits +
                          (sizeof(mphf_configuration) + sizeof(k) + sizeof(m) + sizeof(mm_seed) +
                           sizeof(nkmers) + sizeof(distinct_minimizers)) *
                              8;
    return total_bit_size;
}

mphf::mm_context_t mphf::query(kmer_t kmer, uint64_t minimizer, uint32_t position) const {
    mm_context_t res;
    uint64_t mp_hash = minimizer_order(minimizer);
    std::pair<MinimizerType, std::size_t> dummy = wtree.rank_of(mp_hash);
    MinimizerType mm_type = dummy.first;
    uint64_t mm_type_rank = dummy.second;
    uint64_t sk_size;
    // std::cerr << "mm hash = " << mp_hash << "\n";
    switch (mm_type) {
        case LEFT:
            // locpres_hash = 0;  // because in the elias-fano global vector left positions are the
            // left-most block starting at the beginning locpres_hash +=
            // sizes_and_positions.access(mm_type_rank);  // number of left-KMERS before our bucket
            // locpres_hash += position;  // add local rank
            res.global_rank =
                sizes_and_positions.access(mm_type_rank);  // number of left-KMERS before our bucket
            res.local_rank = position;
            res.type = LEFT;
            // std::cerr << "[LEFT] rank = " << mm_type_rank << ", ";
            // std::cerr << "global shift = " << res.global_rank << ", local shift = " <<
            // res.local_rank;
            break;
        case RIGHT_OR_COLLISION:
            sk_size = sizes_and_positions.diff(right_coll_sizes_start + mm_type_rank);
            if (sk_size == 0) {
                res.global_rank = sizes_and_positions.access(
                    none_pos_start);  // prefix sum of all sizes (sizes of collisions are 0)
                res.local_rank = fallback_kmer_order(kmer);
                res.type = NONE + 1;
                // std::cerr << "[COLLISION] rank = " << none_pos_start << ", global shift = " <<
                // res.global_rank << ", "; std::cerr << "local shift = " << sk_size;
            } else {
                // locpres_hash = sizes_and_positions.access(right_coll_sizes_start + mm_type_rank);
                // // global shift locpres_hash += k - m - p;  // local shift
                res.global_rank = sizes_and_positions.access(right_coll_sizes_start +
                                                             mm_type_rank);  // global shift
                res.local_rank = k - m - position;                           // local shift
                res.type = RIGHT_OR_COLLISION;  // in this case it is only RIGHT
                // std::cerr << "[RIGHT] rank = " << right_coll_sizes_start + mm_type_rank << ", ";
                // std::cerr << "global shift = " << res.global_rank << ", ";
                // std::cerr << "local shift = " << res.local_rank;
            }
            break;
        case MAXIMAL:  // easy case
            // locpres_hash = (k - m + 1) * mm_type_rank + position;  // all maximal k-mer hashes
            // are < than those of all the other types
            res.global_rank = (k - m + 1) * mm_type_rank;
            res.local_rank = position;
            res.type = MAXIMAL;
            // std::cerr << "[MAXIMAL] rank = " << mm_type_rank << ", ";
            // std::cerr << "global shift = " << res.global_rank << ", local shift = " <<
            // res.local_rank;
            break;
        case NONE:
            // locpres_hash = sizes_and_positions.access(none_sizes_start + mm_type_rank);  //
            // prefix sum of sizes locpres_hash += sk_size - position;  // position in the first
            // k-mer - actual position = local shift
            res.global_rank =
                sizes_and_positions.access(none_sizes_start + mm_type_rank);  // prefix sum of sizes
            sk_size = sizes_and_positions.diff(none_pos_start + mm_type_rank);  // p1 actually
            res.local_rank = sk_size - position;
            res.type = NONE;
            // std::cerr << "[NONE] rank = " << none_sizes_start + mm_type_rank << " = " <<
            // none_sizes_start << " + " << mm_type_rank << ", "; std::cerr << "global rank = " <<
            // res.global_rank << ", "; std::cerr << "local shift = " << res.local_rank << ", p1 = "
            // << sk_size << ", p = " << position;
            break;
        default:
            throw std::runtime_error("Unrecognized minimizer type");
    }
    if (mm_type != MAXIMAL) {
        // locpres_hash += (k - m + 1) * n_maximal;  // shift of the maximal k-mers
        res.global_rank += (k - m + 1) * n_maximal;  // shift of the maximal k-mers
    }
    res.hval = res.global_rank + res.local_rank;
    // if (false) {
    //     std::cerr << ", ";
    //     std::string explicit_kmer(contig, i, k);
    //     std::cerr << explicit_kmer << ", ";
    //     std::cerr << kmer << ", ";
    //     std::cerr << "minimizer = " << mm << ", ";
    //     std::cerr << "mm pos = " << p << ", ";
    //     std::cerr << "hash = " << locpres_hash << "\n";
    // }
    return res;
}

void mphf::print_statistics() const noexcept {
    auto mm_mphf_size_bits = minimizer_order.num_bits();
    auto triplet_tree_size_bits = wtree.num_bits();
    auto elias_sequence_size_bits = (sizeof(n_maximal) + sizeof(right_coll_sizes_start) +
                                     sizeof(none_sizes_start) + sizeof(none_pos_start)) *
                                        8 +
                                    sizes_and_positions.num_bits();
    auto kmer_mphf_size_bits = fallback_kmer_order.num_bits();
    auto total_bit_size = mm_mphf_size_bits + triplet_tree_size_bits + elias_sequence_size_bits +
                          kmer_mphf_size_bits +
                          (sizeof(mphf_configuration) + sizeof(k) + sizeof(m) + sizeof(mm_seed) +
                           sizeof(nkmers) + sizeof(distinct_minimizers)) *
                              8;
    std::cerr << "Total number of k-mers: " << nkmers << "\n";
    std::cerr << "Minimizer MPHF size in bits : " << mm_mphf_size_bits << " ("
              << static_cast<double>(mm_mphf_size_bits) / total_bit_size * 100 << "%)\n";
    std::cerr << "\t = " << static_cast<double>(mm_mphf_size_bits) / minimizer_order.num_keys()
              << " bits/minimizer\n\n";
    std::cerr << "Wavelet tree size in bits : " << triplet_tree_size_bits << " ("
              << static_cast<double>(triplet_tree_size_bits) / total_bit_size * 100 << "%)\n";
    std::cerr << "\t = " << static_cast<double>(triplet_tree_size_bits) / minimizer_order.num_keys()
              << " bits/minimizer\n\n";
    std::cerr << "Compressed arrays (EF) : " << elias_sequence_size_bits << " ("
              << static_cast<double>(elias_sequence_size_bits) / total_bit_size * 100 << "%)\n";
    std::cerr << "\t = "
              << static_cast<double>(elias_sequence_size_bits) / sizes_and_positions.size()
              << " bits/offset\n\n";
    std::cerr << "Fallback MPHF : " << kmer_mphf_size_bits << " ("
              << static_cast<double>(kmer_mphf_size_bits) / total_bit_size * 100 << "%)\n";
    std::cerr << "\t = "
              << static_cast<double>(kmer_mphf_size_bits) / fallback_kmer_order.num_keys()
              << " bits/kmer\n\n";
    std::cerr << "Total size in bits : " << total_bit_size << "\n";
    std::cerr << "\tequivalent to : " << static_cast<double>(total_bit_size) / nkmers
              << " bits/k-mer\n";
    std::cerr << "\n";
}

std::ostream& operator<<(std::ostream& out, mphf const& hf) {
    out << "k = " << static_cast<uint32_t>(hf.k) << "\n";
    out << "m = " << static_cast<uint32_t>(hf.m) << "\n";
    out << "minimizer seed = " << hf.mm_seed << "\n";
    out << "number of k-mers = " << hf.nkmers << "\n";
    out << "distinct minimizers = " << hf.distinct_minimizers << "\n";
    out << "maximal super-k-mers = " << hf.n_maximal << "\n";
    out << "starting index of right and collision sizes = " << hf.right_coll_sizes_start << "\n";
    out << "starting index of none sizes = " << hf.none_sizes_start << "\n";
    out << "starting index of none positions = " << hf.none_pos_start << "\n";
    return out;
}

mphf_alt::mphf_alt() : k(0), m(0), mm_seed(0), nkmers(0), distinct_minimizers(0) {
    mphf_configuration.minimal_output = true;
    mphf_configuration.seed = constants::seed;
    mphf_configuration.c = constants::c;
    mphf_configuration.alpha = 0.94;
    mphf_configuration.verbose_output = false;
    mphf_configuration.num_threads = 0;
    mphf_configuration.tmp_dir = "";
};

mphf_alt::mphf_alt(uint8_t klen, uint8_t mm_size, uint64_t seed, uint64_t total_number_of_kmers,
                   uint8_t nthreads, std::string temporary_directory, bool verbose)
    : k(klen), m(mm_size), mm_seed(seed), nkmers(total_number_of_kmers), distinct_minimizers(0) {
    mphf_configuration.minimal_output = true;
    mphf_configuration.seed = constants::seed;
    mphf_configuration.c = constants::c;
    mphf_configuration.alpha = 0.94;
    mphf_configuration.verbose_output = verbose;
    mphf_configuration.num_threads = nthreads;
    if (temporary_directory != "") {
        mphf_configuration.tmp_dir = temporary_directory;
        essentials::create_directory(temporary_directory);
    }
}

std::vector<uint64_t> mphf_alt::build_index(std::vector<mm_triplet_t>& minimizers) {
    std::sort(minimizers.begin(), minimizers.end());
    distinct_minimizers = 0;
    for (auto it = minimizers.begin(), prev = minimizers.begin(); it != minimizers.end(); ++it) {
        // FIXME get number of distinct minimizer -> rework pthash interface for this
        if (prev->itself != it->itself) {
            ++distinct_minimizers;
            prev = it;
        }
    }
    if (minimizers.size()) ++distinct_minimizers;
    auto begin = vector_mm_triplet_to_pthash_itr_adapter(minimizers.begin(), minimizers.end());
    minimizer_order.build_in_external_memory(begin, distinct_minimizers, mphf_configuration);

    auto mphf_compare = [this](mm_triplet_t const& a, mm_triplet_t const& b) {
        return minimizer_order(a.itself) < minimizer_order(b.itself);
    };
    std::sort(minimizers.begin(), minimizers.end(), mphf_compare);

    std::vector<uint64_t> colliding_minimizers;
    std::vector<uint32_t> pos;
    std::vector<uint32_t> siz;
    for (std::size_t i = 0; i < minimizers.size();) {
        std::size_t j;
        for (j = i; j < minimizers.size() && minimizers[i].itself == minimizers[j].itself; ++j) {}
        if (j > (i + 1)) {
            colliding_minimizers.push_back(minimizers[i].itself);
            pos.push_back(0);
            siz.push_back(0);

            // for (std::size_t l = i; l < j; ++l) {
            //     minimizers[l].p1 = 0;
            //     minimizers[l].size = 0;
            // }
        } else {
            // for (; i < j; ++i) {
            pos.push_back(minimizers[i].p1);
            siz.push_back(minimizers[i].size);
        }
        i = j;
    }

    positions.encode(pos.begin(), pos.size());
    sizes.encode(siz.begin(), siz.size());
    num_kmers_in_main_index = sizes.access(sizes.size() - 1);
    return colliding_minimizers;
}

void mphf_alt::build_fallback_mphf(std::vector<kmer_t>& colliding_kmers) {
#ifndef NDEBUG
    std::sort(colliding_kmers.begin(), colliding_kmers.end());
    auto it = std::unique(colliding_kmers.begin(), colliding_kmers.end());
    assert(it !=
           colliding_kmers.end());  // if false: there are some duplicates in unbucketable_kmers
#endif
    fallback_kmer_order.build_in_external_memory(colliding_kmers.begin(), colliding_kmers.size(),
                                                 mphf_configuration);
}

uint64_t mphf_alt::get_minimizer_L0() const noexcept { return distinct_minimizers; }

uint64_t mphf_alt::get_kmer_count() const noexcept { return nkmers; }

uint64_t mphf_alt::num_bits() const noexcept {
    auto mm_mphf_size_bits = minimizer_order.num_bits();
    auto positions_size_bits = positions.num_bits();
    auto sizes_size_bits = sizes.num_bits();
    auto kmer_mphf_size_bits = fallback_kmer_order.num_bits();
    auto total_bit_size = mm_mphf_size_bits + positions_size_bits + sizes_size_bits +
                          kmer_mphf_size_bits +
                          (sizeof(mphf_configuration) + sizeof(k) + sizeof(m) + sizeof(mm_seed) +
                           sizeof(nkmers) + sizeof(distinct_minimizers)) *
                              8;
    return total_bit_size;
}

mphf_alt::mm_context_t mphf_alt::query(kmer_t kmer, uint64_t minimizer, uint32_t position) const {
    mm_context_t res;
    uint64_t index = minimizer_order(minimizer);
    auto [val1, val2] = sizes.pair(index);
    uint64_t size = val2 - val1;
    if (size == 0) {
        res.hval = num_kmers_in_main_index + fallback_kmer_order(kmer);
        res.collision = true;
        return res;
    }
    auto p1 = positions.diff(index);
    res.hval = val1 + p1 - position;
    res.collision = false;
    return res;
}

void mphf_alt::print_statistics() const noexcept {
    auto mm_mphf_size_bits = minimizer_order.num_bits();
    auto positions_size_bits = positions.num_bits();
    auto sizes_size_bits = sizes.num_bits();
    auto kmer_mphf_size_bits = fallback_kmer_order.num_bits();
    auto total_bit_size = mm_mphf_size_bits + positions_size_bits + sizes_size_bits +
                          kmer_mphf_size_bits +
                          (sizeof(mphf_configuration) + sizeof(k) + sizeof(m) + sizeof(mm_seed) +
                           sizeof(nkmers) + sizeof(distinct_minimizers)) *
                              8;
    std::cerr << "Total number of k-mers: " << nkmers << "\n";
    std::cerr << "Minimizer MPHF size in bits : " << mm_mphf_size_bits << " ("
              << static_cast<double>(mm_mphf_size_bits) / total_bit_size * 100 << "%)\n";
    std::cerr << "\t = " << static_cast<double>(mm_mphf_size_bits) / minimizer_order.num_keys()
              << " bits/minimizer\n\n";
    std::cerr << "positions EF sequence size in bits : " << positions_size_bits << " ("
              << static_cast<double>(positions_size_bits) / total_bit_size * 100 << "%)\n";
    std::cerr << "\t = " << static_cast<double>(positions_size_bits) / minimizer_order.num_keys()
              << " bits/minimizer\n\n";
    std::cerr << "sizes EF sequence size in bits : " << sizes_size_bits << " ("
              << static_cast<double>(sizes_size_bits) / total_bit_size * 100 << "%)\n";
    std::cerr << "Fallback MPHF : " << kmer_mphf_size_bits << " ("
              << static_cast<double>(kmer_mphf_size_bits) / total_bit_size * 100 << "%)\n";
    std::cerr << "\t = "
              << static_cast<double>(kmer_mphf_size_bits) / fallback_kmer_order.num_keys()
              << " bits/kmer\n\n";
    std::cerr << "Total size in bits : " << total_bit_size << "\n";
    std::cerr << "\tequivalent to : " << static_cast<double>(total_bit_size) / nkmers
              << " bits/k-mer\n";
    std::cerr << "\n";
}

std::ostream& operator<<(std::ostream& out, mphf_alt const& hf) {
    out << "k = " << static_cast<uint32_t>(hf.k) << "\n";
    out << "m = " << static_cast<uint32_t>(hf.m) << "\n";
    out << "minimizer seed = " << hf.mm_seed << "\n";
    out << "number of k-mers = " << hf.nkmers << "\n";
    out << "distinct minimizers = " << hf.distinct_minimizers << "\n";
    return out;
}

}  // namespace lphash