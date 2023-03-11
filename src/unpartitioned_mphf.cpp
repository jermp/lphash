#include <zlib.h>
extern "C" {
#include "../external/kseq.h"
}
#include "../include/unpartitioned_mphf.hpp"
#include "../include/minimizer.hpp"

KSEQ_INIT(gzFile, gzread)

namespace lphash {

mphf_alt::mphf_alt() : k(0), m(0), mm_seed(0), nkmers(0), distinct_minimizers(0), max_ram(false) {
    mphf_configuration.minimal_output = true;
    mphf_configuration.seed = constants::default_pthash_seed;
    mphf_configuration.c = constants::c;
    mphf_configuration.alpha = 0.94;
    mphf_configuration.verbose_output = false;
    mphf_configuration.num_threads = 0;
    mphf_configuration.tmp_dir = "";
    mphf_configuration.ram = 8 * essentials::GB;
};

void mphf_alt::build(configuration const& config, std::ostream& res_strm) {
    constexpr bool canonical = false;
    k = config.k;
    m = config.m;
    mm_seed = config.mm_seed;
    nkmers = 0;
    distinct_minimizers = 0;
    max_ram = config.max_memory;
    mphf_configuration.minimal_output = true;
    mphf_configuration.seed = constants::default_pthash_seed;
    mphf_configuration.c = config.c;
    mphf_configuration.alpha = 0.94;
    mphf_configuration.verbose_output = config.verbose;
    mphf_configuration.num_threads = config.num_threads;
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
        mphf_configuration.tmp_dir, get_group_id());
    if ((fp = gzopen(config.input_filename.c_str(), "r")) == NULL)
        throw std::runtime_error("Unable to open the input file " + config.input_filename + "\n");
    uint64_t id = 0;
    seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        auto n = minimizer::from_string<pthash::murmurhash2_64>(
            seq->seq.s, seq->seq.l, k, m, mm_seed, canonical, id,
            all_minimizers);  // non-canonical minimizers for now
        nkmers += n;
        ++total_contigs;
        check_total_kmers += seq->seq.l - k + 1;
    }
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
        uint64_t pos_sum, size_sum;
        pos_sum = size_sum = 0;
        for (auto itr = unique_mms.cbegin(); itr != unique_mms.cend(); ++itr) {
            auto triplet = *itr;
            triplet.itself = get_minimizer_order(triplet.itself);
            mm_sorted_by_mphf.push_back(triplet);
            pos_sum += triplet.p1;
            size_sum += triplet.size;
        }
        explicit_garbage_collect(std::move(unique_mms));
        auto itr = mm_sorted_by_mphf.cbegin();
        build_pos_index(itr, mm_sorted_by_mphf.size(), pos_sum);
        itr = mm_sorted_by_mphf.cbegin();
        build_size_index(itr, mm_sorted_by_mphf.size(), size_sum);
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
            throw std::runtime_error("Unable to open the input file " + config.input_filename +
                                     " a seconf time\n");
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
    if (total_contigs >= 1) --total_contigs;
    res_strm << config.input_filename << "," << static_cast<uint32_t>(k) << ","
             << static_cast<uint32_t>(m) << ","
             << static_cast<double>(total_colliding_minimizers) / distinct_minimizers << ","
             << 2.0 / ((k - m + 1) + 1) << "," << static_cast<double>(total_minimizers) / nkmers
             << "," << static_cast<double>(total_contigs) / nkmers << ","
             << static_cast<double>(num_bits()) / get_kmer_count();
    res_strm << "\n";
}

void mphf_alt::build_minimizers_mphf(
    external_memory_vector<mm_triplet_t, false>::const_iterator& mm_itr,
    std::size_t number_of_distinct_minimizers) {
    mm_itr_t dummy_itr(mm_itr);
    distinct_minimizers = number_of_distinct_minimizers;
    minimizer_order.build_in_external_memory(dummy_itr, distinct_minimizers, mphf_configuration);
}

void mphf_alt::build_fallback_mphf(external_memory_vector<kmer_t, false>::const_iterator& km_itr,
                                   std::size_t number_of_colliding_kmers) {
    km_itr_t dummy_itr(km_itr);
    fallback_kmer_order.build_in_external_memory(dummy_itr, number_of_colliding_kmers,
                                                 mphf_configuration);
}

void mphf_alt::build_pos_index(external_memory_vector<mm_triplet_t>::const_iterator& mm_itr,
                               [[maybe_unused]] std::size_t number_of_distinct_minimizers,
                               uint64_t pos_sum) {
    assert(distinct_minimizers == number_of_distinct_minimizers);
    pos_itr_t dummy_itr(mm_itr);
    cumulative_iterator c_itr(dummy_itr);
    positions.encode(c_itr, distinct_minimizers, pos_sum);
}

void mphf_alt::build_size_index(external_memory_vector<mm_triplet_t>::const_iterator& mm_itr,
                                [[maybe_unused]] std::size_t number_of_distinct_minimizers,
                                uint64_t size_sum) {
    assert(distinct_minimizers == number_of_distinct_minimizers);
    size_itr_t dummy_itr(mm_itr);
    cumulative_iterator c_itr(dummy_itr);
    sizes.encode(c_itr, distinct_minimizers, size_sum);
    num_kmers_in_main_index = sizes.access(sizes.size() - 1);
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

uint64_t mphf_alt::get_minimizer_order(uint64_t mm) const { return minimizer_order(mm); }

mphf_alt::mm_context_t mphf_alt::query(kmer_t kmer, uint64_t minimizer, uint32_t position) const {
    mm_context_t res;
    uint64_t index = minimizer_order(minimizer);
    auto [val1, val2] = sizes.pair(index);
    assert(val2 >= val1);
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
    std::cerr << "Total number of k-mers belonging to ambiguous minimizers: "
              << fallback_kmer_order.num_keys() << "\n";
    std::cerr << "xi = " << static_cast<double>(fallback_kmer_order.num_keys()) / nkmers << "\n";
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