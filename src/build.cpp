extern "C" {
#include "../external/kseq.h"
}
#include <zlib.h>
#include "../include/constants.hpp"
#include "../include/mphf.hpp"
#include "parser_config.hpp"
#include "minimizer.hpp"

using namespace lphash;

KSEQ_INIT(gzFile, gzread)

int main(int argc, char* argv[]) {
    cmd_line_parser::parser parser(argc, argv);
    try {
        parser = get_build_parser(argc, argv);
    } catch (const ParseError& e) { return 1; }

    std::string input_filename = parser.get<std::string>("input_filename");

    uint64_t k = parser.get<uint64_t>("k");
    if (k > constants::max_k) {
        std::cerr << "k cannot be larger than " + std::to_string(constants::max_k) + "\n";
        return 1;
    }

    uint64_t m = parser.get<uint64_t>("m");
    if (m > k) {
        std::cerr << "m cannot be larger than k\n";
        return 1;
    }

    uint64_t mm_seed = constants::default_seed;
    if (parser.parsed("seed")) mm_seed = parser.get<uint64_t>("seed");

    uint64_t num_threads = constants::default_num_threads;
    if (parser.parsed("threads")) num_threads = parser.get<uint64_t>("threads");

    std::string tmp_dirname = constants::default_tmp_dirname;
    if (parser.parsed("tmp_dirname")) {
        tmp_dirname = parser.get<std::string>("tmp_dirname");
        essentials::create_directory(tmp_dirname);
    }

    double c = constants::c;
    if (parser.parsed("c")) {
        c = parser.get<double>("c");
        if (c > 10.0 || c < 3.0) {
            std::cerr << "3.0 <= c <= 10.0\n";
            return 1;
        }
    }

    uint64_t max_memory = 8;
    if (parser.parsed("max-memory")) {
        max_memory = parser.get<uint64_t>("max-memory");
        if (max_memory > 255) {
            std::cerr << "The maximum allowed amount of ram is 255GB\n" << std::endl;
            return 1;
        }
    }

    bool check = false;
    bool verbose = false;
    if (parser.parsed("check")) check = parser.get<bool>("check");
    if (parser.parsed("verbose")) verbose = parser.get<bool>("verbose");

    uint64_t total_kmers = 0;
    uint64_t check_total_kmers = 0;
    uint64_t total_contigs = 0;
    std::cerr << "Part 1: file reading and info gathering\n";
    sorted_external_vector<mm_record_t> all_minimizers(
        uint64_t(max_memory) * essentials::GB,
        [](mm_record_t const& a, mm_record_t const& b) { return a.itself < b.itself; }, tmp_dirname,
        get_group_id());

    gzFile fp = NULL;
    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
        std::cerr << "Unable to open the input file " << input_filename << "\n";
        return 2;
    }

    constexpr bool canonical = false;
    kseq_t* seq = kseq_init(fp);
    uint64_t id = 0;
    while (kseq_read(seq) >= 0) {
        uint64_t n =
            minimizer::from_string<hash64>(seq->seq.s, seq->seq.l, k, m, mm_seed, canonical, id,
                                           all_minimizers);  // non-canonical minimizers for now
        total_kmers += n;
        ++total_contigs;
        check_total_kmers += seq->seq.l - k + 1;
    }
    if (total_contigs > 0) --total_contigs;
    if (seq) kseq_destroy(seq);
    gzclose(fp);
    uint64_t total_minimizers = all_minimizers.size();
    assert(total_kmers == check_total_kmers);

    std::cerr << "Part 2: build MPHF\n";
    auto [unique_mms, coll_ids] =
        minimizer::classify(std::move(all_minimizers), max_memory, tmp_dirname);
    mphf f(k, m, mm_seed, total_kmers, c, num_threads, max_memory, tmp_dirname, verbose);
    {
        auto itr = unique_mms.cbegin();
        f.build_minimizers_mphf(itr, unique_mms.size());
    }

    std::cerr << "Part 3: build inverted index\n";
    {
        sorted_external_vector<mm_triplet_t> mm_sorted_by_mphf(
            uint64_t(max_memory) * essentials::GB,
            [](mm_triplet_t const& a, mm_triplet_t const& b) { return a.itself < b.itself; },
            tmp_dirname, get_group_id());
        for (auto itr = unique_mms.cbegin(); itr != unique_mms.cend(); ++itr) {
            auto triplet = *itr;
            triplet.itself = f.get_minimizer_order(triplet.itself);
            mm_sorted_by_mphf.push_back(triplet);
        }
        explicit_garbage_collect(std::move(unique_mms));
        auto itr = mm_sorted_by_mphf.cbegin();
        f.build_inverted_index(itr, mm_sorted_by_mphf.size());
    }
    uint64_t total_distinct_minimizers = f.get_minimizer_L0();
    uint64_t total_colliding_minimizers = coll_ids.size();

    std::cerr << "Part 4: build fallback MPHF\n";
    {  // garbage collector for unbucketable_kmers
        sorted_external_vector<kmer_t> unbucketable_kmers(
            uint64_t(max_memory) * essentials::GB,
            []([[maybe_unused]] kmer_t const& a, [[maybe_unused]] kmer_t const& b) {
                return false;
            },
            tmp_dirname, get_group_id());
        std::unordered_map<uint64_t, uint64_t> stats;
        if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
            std::cerr << "Unable to open the input file a second time" << input_filename << "\n";
            return 2;
        }
        id = 0;
        auto start = coll_ids.cbegin();
        auto stop = coll_ids.cend();
        seq = kseq_init(fp);
        while (kseq_read(seq) >= 0) {
            minimizer::get_colliding_kmers<hash64>(seq->seq.s, seq->seq.l, k, m, mm_seed, canonical,
                                                   start, stop, id, unbucketable_kmers, stats);
        }
        if (seq) kseq_destroy(seq);
        gzclose(fp);
        explicit_garbage_collect(std::move(coll_ids));
        {
            auto itr = unbucketable_kmers.cbegin();
            f.build_fallback_mphf(itr, unbucketable_kmers.size());
        }
    }

    if (parser.parsed("output_filename")) {
        auto output_filename = parser.get<std::string>("output_filename");
        std::cerr << "Saving data structure to disk...  ";
        essentials::save(f, output_filename.c_str());
        std::cerr << "DONE\n";
    }

    if (check) {
        std::cerr << "Checking...\n";
        if (parser.parsed("output_filename")) {
            mphf loaded;
            uint64_t num_bytes_read =
                essentials::load(loaded, parser.get<std::string>("output_filename").c_str());
            std::cerr << "[Info] Loaded " << num_bytes_read * 8 << " bits\n";
            pthash::bit_vector_builder population(
                loaded.get_kmer_count());  // bitvector for checking perfection and minimality
            if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
                std::cerr << "Unable to open the input file a second time" << input_filename
                          << "\n";
                return 2;
            }
            seq = kseq_init(fp);
            while (check && kseq_read(seq) >= 0) {
                std::string contig = std::string(seq->seq.s);
                check = check_collisions(loaded, contig, canonical, population);
                if (check) check = check_streaming_correctness(loaded, contig, canonical);
            }
            if (seq) kseq_destroy(seq);
            gzclose(fp);
            check = check && check_perfection(loaded, population);
        } else {
            pthash::bit_vector_builder population(f.get_kmer_count());
            if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
                std::cerr << "Unable to open the input file a second time" << input_filename
                          << "\n";
                return 2;
            }
            seq = kseq_init(fp);
            while (check && kseq_read(seq) >= 0) {
                std::string contig = std::string(seq->seq.s);
                check = check_collisions(f, contig, canonical, population);
                if (check) check = check_streaming_correctness(f, contig, canonical);
            }
            if (seq) kseq_destroy(seq);
            gzclose(fp);
            check = check && check_perfection(f, population);
        }
    }

    if (verbose) {
        std::cerr << "Statistics:\n";
        f.print_statistics();
    }

    std::cout << input_filename << "," << static_cast<uint32_t>(k) << ","
              << static_cast<uint32_t>(m) << ","
              << static_cast<double>(total_colliding_minimizers) / total_distinct_minimizers << ","
              << 2.0 / ((k - m + 1) + 1) << ","
              << static_cast<double>(total_minimizers) / total_kmers << ","
              << static_cast<double>(total_contigs) / total_kmers << ","
              << static_cast<double>(f.num_bits()) / f.get_kmer_count();
    std::cout << "\n";

    return 0;
}