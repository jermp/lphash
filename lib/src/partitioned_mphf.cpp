#include <zlib.h>
extern "C" {
#include "../../external/kseq.h"
}
#include "../include/partitioned_mphf.hpp"
#include "../include/minimizer.hpp"

KSEQ_INIT(gzFile, gzread)

namespace lphash {
namespace mphf {

partitioned::partitioned()
    : interface()
    , n_maximal(0)
    , right_coll_sizes_start(0)
    , none_sizes_start(0)
    , none_pos_start(0) {
    mphf_configuration.minimal_output = true;
    mphf_configuration.seed = 0;
    mphf_configuration.c = 0;
    mphf_configuration.alpha = 0;
    mphf_configuration.verbose_output = false;
    mphf_configuration.num_threads = 0;
    mphf_configuration.tmp_dir = "";
    mphf_configuration.ram = 0;
};

void partitioned::build(configuration const& config, std::ostream& res_strm) {
    constexpr bool canonical = false;
    k = config.k;
    m = config.m;
    mm_seed = config.mm_seed;
    nkmers = 0;
    distinct_minimizers = 0;
    n_maximal = 0;
    right_coll_sizes_start = 0;
    none_sizes_start = 0;
    none_pos_start = 0;
    max_ram = config.max_memory;
    mphf_configuration.minimal_output = true;
    mphf_configuration.seed = config.pt_seed;
    mphf_configuration.c = config.c;
    mphf_configuration.alpha = 0.94;
    mphf_configuration.verbose_output = config.verbose;
    mphf_configuration.num_threads = config.num_threads;
    mphf_configuration.num_partitions = config.num_threads; // FIXME add proper option
    mphf_configuration.ram = static_cast<uint64_t>(max_ram) * essentials::GB;
    mphf_configuration.tmp_dir = config.tmp_dirname;

    uint64_t check_total_kmers = 0;
    uint64_t total_contigs = 0;
    uint64_t total_minimizers = 0;
    uint64_t total_colliding_minimizers = 0;

    gzFile fp = nullptr;
    kseq_t* seq = nullptr;
    if (config.verbose) std::cerr << "Part 1: file reading and info gathering\n";
    external_memory_vector<mm_record_t> all_minimizers(
        uint64_t(max_ram) * essentials::GB,
        [](mm_record_t const& a, mm_record_t const& b) { return a.itself < b.itself; },
        config.tmp_dirname, get_group_id());
    if ((fp = gzopen(config.input_filename.c_str(), "r")) == NULL)
        throw std::runtime_error("Unable to open the input file " + config.input_filename + "\n");
    seq = kseq_init(fp);
    uint64_t id = 0;
    while (kseq_read(seq) >= 0) {
        uint64_t n = minimizer::from_string<pthash::murmurhash2_64>(
            seq->seq.s, seq->seq.l, k, m, mm_seed, canonical, id,
            all_minimizers);  // non-canonical minimizers for now
        nkmers += n;
        ++total_contigs;
        check_total_kmers += seq->seq.l - k + 1;
    }
    if (total_contigs > 0) --total_contigs;
    if (seq) kseq_destroy(seq);
    gzclose(fp);
    total_minimizers = all_minimizers.size();
    assert(nkmers == check_total_kmers);

    if (config.verbose) std::cerr << "Part 2: build MPHF\n";
    auto [unique_mms, coll_ids] =
        minimizer::classify(std::move(all_minimizers), max_ram, mphf_configuration.tmp_dir);
    {
        auto itr = unique_mms.cbegin();
        build_minimizers_mphf(itr, unique_mms.size());
    }

    if (config.verbose) std::cerr << "Part 3: build inverted index\n";
    {
        external_memory_vector<mm_triplet_t> mm_sorted_by_mphf(
            uint64_t(max_ram) * essentials::GB,
            [](mm_triplet_t const& a, mm_triplet_t const& b) { return a.itself < b.itself; },
            mphf_configuration.tmp_dir, get_group_id());
        for (auto itr = unique_mms.cbegin(); itr != unique_mms.cend(); ++itr) {
            auto triplet = *itr;
            triplet.itself = get_minimizer_order(triplet.itself);
            mm_sorted_by_mphf.push_back(triplet);
        }
        explicit_garbage_collect(std::move(unique_mms));
        auto itr = mm_sorted_by_mphf.cbegin();
        build_inverted_index(itr, mm_sorted_by_mphf.size());
    }

    total_colliding_minimizers = coll_ids.size();

    if (config.verbose) std::cerr << "Part 4: build fallback MPHF\n";
    {  // garbage collector for unbucketable_kmers
        external_memory_vector<kmer_t, false> unbucketable_kmers(
            uint64_t(max_ram) * essentials::GB,
            // []([[maybe_unused]] kmer_t const& a, [[maybe_unused]] kmer_t const& b) {
            //     return false;
            // },
            mphf_configuration.tmp_dir, get_group_id());
        if ((fp = gzopen(config.input_filename.c_str(), "r")) == NULL)
            throw std::runtime_error("Unable to open the input file a second time " +
                                     config.input_filename + "\n");
        id = 0;
        auto start = coll_ids.cbegin();
        auto stop = coll_ids.cend();
        seq = kseq_init(fp);
        while (kseq_read(seq) >= 0) {
            minimizer::get_colliding_kmers<pthash::murmurhash2_64>(seq->seq.s, seq->seq.l, k, m,
                                                                   mm_seed, canonical, start, stop,
                                                                   id, unbucketable_kmers);
        }
        if (seq) kseq_destroy(seq);
        gzclose(fp);
        explicit_garbage_collect(std::move(coll_ids));
        {
            auto itr = unbucketable_kmers.cbegin();
            build_fallback_mphf(itr, unbucketable_kmers.size());
        }
    }
    res_strm << config.input_filename << "," << static_cast<uint32_t>(k) << ","
             << static_cast<uint32_t>(m) << ","
             << static_cast<double>(total_colliding_minimizers) / distinct_minimizers << ","
             << 2.0 / ((k - m + 1) + 1) << "," << static_cast<double>(total_minimizers) / nkmers
             << "," << static_cast<double>(total_contigs) / nkmers << ","
             << static_cast<double>(num_bits()) / nkmers;
    res_strm << "\n";
}

// void partitioned::build_minimizers_mphf(
//     external_memory_vector<mm_triplet_t, false>::const_iterator& mm_itr,
//     std::size_t number_of_distinct_minimizers) {
//     mm_itr_t dummy_itr(mm_itr);
//     distinct_minimizers = number_of_distinct_minimizers;
//     minimizer_order.build_in_external_memory(std::move(dummy_itr), distinct_minimizers, mphf_configuration);
// }

// void partitioned::build_fallback_mphf(external_memory_vector<kmer_t, false>::const_iterator& km_itr,
//                                std::size_t number_of_colliding_kmers) {
//     km_itr_t dummy_itr(km_itr);
//     fallback_kmer_order.build_in_external_memory(dummy_itr, number_of_colliding_kmers,
//                                                  mphf_configuration);
// }

void partitioned::build_inverted_index(external_memory_vector<mm_triplet_t>::const_iterator& mm_itr,
                                [[maybe_unused]] std::size_t number_of_distinct_minimizers) {
    assert(distinct_minimizers == number_of_distinct_minimizers);
    n_maximal = 0;
    std::size_t colliding_minimizers = 0;
    uint64_t universe = 0;
    quartet_wtree_builder wtb(distinct_minimizers);

    uint64_t array_mem = max_ram * essentials::GB / 4;
    array_mem = array_mem < 4000000 ? 4000000 : array_mem;
    external_memory_vector<uint64_t, false> left_positions(array_mem, mphf_configuration.tmp_dir,
                                                           get_group_id());
    external_memory_vector<uint64_t, false> right_or_collision_sizes(
        array_mem, mphf_configuration.tmp_dir, get_group_id());
    external_memory_vector<uint64_t, false> none_sizes(array_mem, mphf_configuration.tmp_dir,
                                                       get_group_id());
    external_memory_vector<uint64_t, false> none_positions(array_mem, mphf_configuration.tmp_dir,
                                                           get_group_id());
    typedef external_memory_vector<uint64_t, false>::const_iterator itr_t;

    for (std::size_t i = 0; i < distinct_minimizers; ++i) {
        mm_triplet_t mm = *mm_itr;
        if (mm.size == 0) {
            wtb.push_back(RIGHT_OR_COLLISION);
            right_or_collision_sizes.push_back(0);
            // universe += 0;
            ++colliding_minimizers;
        } else {
            if (mm.p1 == k - m) {
                if (mm.size == k - m + 1) {
                    wtb.push_back(MAXIMAL);
                    ++n_maximal;
                } else {
                    wtb.push_back(RIGHT_OR_COLLISION);
                    right_or_collision_sizes.push_back(mm.size);
                    universe += mm.size;
                }
            } else {
                if (mm.p1 == mm.size - 1) {
                    wtb.push_back(LEFT);
                    left_positions.push_back(
                        mm.p1 + 1);  // +1 because with p1 == 0 we have 1 k-mer in the prefix sum
                    universe += mm.p1 + 1;
                } else {
                    wtb.push_back(NONE);
                    none_positions.push_back(mm.p1);  // here we do not have +1 since p1 != 0 by
                                                      // itself (NONE type of super-k-mer)
                    none_sizes.push_back(mm.size);
                    universe += mm.p1 + mm.size;
                }
            }
        }
        ++mm_itr;
    }
    assert(none_positions.size() == none_sizes.size());

    wtree.build(wtb);

    right_coll_sizes_start = left_positions.size();
    none_sizes_start = right_coll_sizes_start + right_or_collision_sizes.size();
    none_pos_start = none_sizes_start + none_sizes.size();

    if (mphf_configuration.verbose_output) {
        // double maximal = static_cast<double>(n_maximal) / distinct_minimizers * 100;
        // double left = static_cast<double>(left_positions.size()) / distinct_minimizers * 100;
        // double right = static_cast<double>(right_or_collision_sizes.size() -
        // colliding_minimizers) /
        //                distinct_minimizers * 100;
        // double none = static_cast<double>(none_positions.size()) / distinct_minimizers * 100;
        // double ambiguous = static_cast<double>(colliding_minimizers) / distinct_minimizers * 100;
        // std::cerr << "Percentage of maximal super-k-mers: " << maximal << "%\n";
        // std::cerr << "Percentage of left-maximal super-k-mers: " << left << "%\n";
        // std::cerr << "Percentage of right-maximal super-k-mers : " << right << "%\n";
        // std::cerr << "Percentage of unclassified super-k-mers: " << none << "%\n";
        // std::cerr << "Percentage of ambiguous minimizers: " << ambiguous << "%\n";
        uint64_t n = distinct_minimizers - colliding_minimizers;
        double maximal = static_cast<double>(n_maximal) / n * 100;
        double left = static_cast<double>(left_positions.size()) / n * 100;
        double right =
            static_cast<double>(right_or_collision_sizes.size() - colliding_minimizers) / n * 100;
        double none = static_cast<double>(none_positions.size()) / n * 100;
        double ambiguous = static_cast<double>(colliding_minimizers) / distinct_minimizers * 100;
        std::cerr << "Percentage of maximal super-k-mers: " << maximal << "%\n";
        std::cerr << "Percentage of left-maximal super-k-mers: " << left << "%\n";
        std::cerr << "Percentage of right-maximal super-k-mers : " << right << "%\n";
        std::cerr << "Percentage of unclassified super-k-mers: " << none << "%\n";
        std::cerr << "Total = " << (maximal + left + right + none) << "%\n";
        std::cerr << "Percentage of ambiguous minimizers: " << ambiguous << "%\n";
    }

    std::vector<std::pair<itr_t, itr_t>> sp;
    // The following is needed to build the pairs in-place, since const_iterator is non
    // copy-constructible
    sp.emplace_back(std::piecewise_construct, std::make_tuple(left_positions.cbegin()),
                    std::make_tuple(left_positions.cend()));
    sp.emplace_back(std::piecewise_construct, std::make_tuple(right_or_collision_sizes.cbegin()),
                    std::make_tuple(right_or_collision_sizes.cend()));
    sp.emplace_back(std::piecewise_construct, std::make_tuple(none_sizes.cbegin()),
                    std::make_tuple(none_sizes.cend()));
    sp.emplace_back(std::piecewise_construct, std::make_tuple(none_positions.cbegin()),
                    std::make_tuple(none_positions.cend()));
    append_iterator sp_itr(sp);
    cumulative_iterator c_itr(sp_itr);
    assert(distinct_minimizers == none_pos_start + n_maximal);
    sizes_and_positions.encode(c_itr, none_pos_start + none_positions.size(), universe);
}

uint64_t partitioned::num_bits() const noexcept {
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

partitioned::mm_context_t partitioned::query(kmer_t kmer, uint64_t minimizer, uint32_t position) const {
    mm_context_t res;
    uint64_t mp_hash = minimizer_order(minimizer);
    auto [mm_type, mm_type_rank] = wtree.rank_of(mp_hash);

    switch (mm_type) {
        case LEFT:
            res.global_rank = sizes_and_positions.access(mm_type_rank) +
                              (k - m + 1) * n_maximal;  // number of left-KMERS before our bucket
            res.local_rank = position;
            res.type = LEFT;
            break;
        case RIGHT_OR_COLLISION: {
            auto [val1, val2] = sizes_and_positions.pair(right_coll_sizes_start + mm_type_rank);
            assert(val2 >= val1);
            uint64_t sk_size = val2 - val1;
            if (sk_size == 0) {
                res.global_rank =
                    sizes_and_positions.access(none_pos_start) +
                    (k - m + 1) * n_maximal;  // prefix sum of all sizes (sizes of collisions are 0)
                res.local_rank = fallback_kmer_order(kmer);
                res.type = NONE + 1;
            } else {
                res.global_rank = val1 + (k - m + 1) * n_maximal;  // global shift
                res.local_rank = k - m - position;                 // local shift
                res.type = RIGHT_OR_COLLISION;                     // in this case it is only RIGHT
            }
        } break;
        case MAXIMAL:  // easy case
            // maximal k-mer hashes are smaller than those of all the other types
            res.global_rank = (k - m + 1) * mm_type_rank;
            res.local_rank = position;
            res.type = MAXIMAL;
            break;
        case NONE: {
            res.global_rank = sizes_and_positions.access(none_sizes_start + mm_type_rank) +
                              (k - m + 1) * n_maximal;
            uint64_t sk_size =
                sizes_and_positions.diff(none_pos_start + mm_type_rank);  // p1 actually
            res.local_rank = sk_size - position;
            res.type = NONE;
        } break;
        default:
            throw std::runtime_error("Unrecognized minimizer type");
    }
    res.hval = res.global_rank + res.local_rank;
    return res;
}

void partitioned::print_statistics() const noexcept {
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
    std::cerr << "Total number of k-mers belonging to ambiguous minimizers: "
              << fallback_kmer_order.num_keys() << "\n";
    std::cerr << "xi = " << static_cast<double>(fallback_kmer_order.num_keys()) / nkmers << "\n";
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

std::ostream& operator<<(std::ostream& out, partitioned const& hf) {
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

} // namespace mphf
}  // namespace lphash