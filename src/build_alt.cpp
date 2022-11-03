extern "C" {
#include "../include/kseq.h"
}
#include <zlib.h>
#include "../include/constants.hpp"
#include "../include/parser_build.hpp"
#include "../include/mphf.hpp"
#include "minimizer.hpp"

// #include "../include/prettyprint.hpp"

using namespace lphash;

KSEQ_INIT(gzFile, gzread)

int main(int argc, char* argv[]) {
    gzFile fp;
    kseq_t* seq;
    uint8_t k, m, nthreads, max_memory;
    uint64_t mm_seed, id;
    std::string tmp_dirname;
    double c;
    bool canonical;
    bool check;
    bool verbose;
    cmd_line_parser::parser parser(argc, argv);
    std::size_t total_kmers, check_total_kmers, total_minimizers, total_contigs, total_colliding_minimizers, total_distinct_minimizers;

    try {
        parser = get_build_parser(argc, argv);
    } catch (const ParseError& e) {
        return 1;
    }
    std::string input_filename = parser.get<std::string>("input_filename");
    k = static_cast<uint8_t>(parser.get<uint32_t>("k"));
    m = static_cast<uint8_t>(parser.get<uint32_t>("m"));
    if (parser.parsed("seed")) mm_seed = parser.get<uint64_t>("seed");
    else mm_seed = 42;
    if (parser.parsed("threads")) nthreads = parser.get<uint32_t>("threads");
    else nthreads = 1;
    if (parser.parsed("tmp_dirname")) tmp_dirname = parser.get<std::string>("tmp_dirname");
    else tmp_dirname = "";
    if (parser.parsed("c")) c = parser.get<double>("c");
    else c = constants::c;
    if (parser.parsed("canonical_parsing")) canonical = parser.get<bool>("canonical_parsing");
    else canonical = false;
    if (parser.parsed("max-memory")) {
        auto val = parser.get<uint64_t>("max-memory");
        if (val > 255) {
            std::cerr << "The maximum allowed amount of ram is 255GB\n" << std::endl;
            return 1;
        }
        max_memory = static_cast<uint8_t>(val);
    } else max_memory = 8;
    if (parser.parsed("check")) check = parser.get<bool>("check");
    else check = false;
    if (parser.parsed("verbose")) verbose = parser.get<bool>("verbose");
    else verbose = false;
    

    if (k > 64) {
        std::cerr << "k cannot be larger than " + std::to_string(k) + "\n";
        return 1;
    }
    if (m > k) {
        std::cerr << "m cannot be larger than k\n";
        return 1;
    }
    if (c > 10 || c < 3) {
        std::cerr << "3 <= c <= 10\n";
        return 1;
    }

    total_kmers = 0;
    total_contigs = 0;
    check_total_kmers = 0;
    std::cerr << "Part 1: file reading and info gathering\n";
    sorted_external_vector<mm_record_t> all_minimizers(uint64_t(max_memory) * essentials::GB, [](mm_record_t const& a, mm_record_t const& b) {return a.itself < b.itself;}, tmp_dirname, get_group_id());
    fp = NULL;
    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
        std::cerr << "Unable to open the input file " << input_filename << "\n";
        return 2;
    }
    id = 0;
    seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        auto n = minimizer::from_string<hash64>(seq->seq.s, seq->seq.l, k, m, mm_seed, canonical, id, all_minimizers);  // non-canonical minimizers for now
        total_kmers += n;
        ++total_contigs;
        check_total_kmers += seq->seq.l - k + 1;
    }
    if (seq) kseq_destroy(seq);
    gzclose(fp);
    total_minimizers = all_minimizers.size();
    assert(total_kmers == check_total_kmers);

    std::cerr << "Part 2: build MPHF\n";
    auto [unique_mms, coll_ids] = minimizer::classify(std::move(all_minimizers), max_memory, tmp_dirname);
    mphf_alt locpres_mphf(k, m, mm_seed, total_kmers, c, nthreads, max_memory, tmp_dirname, verbose);
    {
        auto itr = unique_mms.cbegin();
        locpres_mphf.build_minimizers_mphf(itr, unique_mms.size());
    }

    std::cerr << "Part 3: build inverted index\n";
    {
        sorted_external_vector<mm_triplet_t> mm_sorted_by_mphf(uint64_t(max_memory) * essentials::GB, [](mm_triplet_t const& a, mm_triplet_t const& b) {return a.itself < b.itself;}, tmp_dirname, get_group_id());
        uint64_t pos_sum, size_sum;
        pos_sum = size_sum = 0;
        for(auto itr = unique_mms.cbegin(); itr != unique_mms.cend(); ++itr) {
            auto triplet = *itr;
            triplet.itself = locpres_mphf.get_minimizer_order(triplet.itself);
            mm_sorted_by_mphf.push_back(triplet);
            pos_sum += triplet.p1;
            size_sum += triplet.size;
        }
        explicit_garbage_collect(std::move(unique_mms));
        auto itr = mm_sorted_by_mphf.cbegin();
        locpres_mphf.build_pos_index(itr, mm_sorted_by_mphf.size(), pos_sum);
        itr = mm_sorted_by_mphf.cbegin();
        locpres_mphf.build_size_index(itr, mm_sorted_by_mphf.size(), size_sum);
    }
    total_distinct_minimizers = locpres_mphf.get_minimizer_L0();
    total_colliding_minimizers = coll_ids.size();

    std::cerr << "Part 4: build fallback MPHF\n";
    {// garbage collector for unbucketable_kmers
        sorted_external_vector<kmer_t> unbucketable_kmers(uint64_t(max_memory) * essentials::GB, []([[maybe_unused]] kmer_t const& a, [[maybe_unused]] kmer_t const& b) {return false;}, tmp_dirname, get_group_id());
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
            minimizer::get_colliding_kmers<hash64>(seq->seq.s, seq->seq.l, k, m, mm_seed, canonical, start, stop, id, unbucketable_kmers, stats);
        }
        if (seq) kseq_destroy(seq);
        gzclose(fp);
        explicit_garbage_collect(std::move(coll_ids));
        {
            auto itr = unbucketable_kmers.cbegin();
            locpres_mphf.build_fallback_mphf(itr, unbucketable_kmers.size());
        }
    }

    if (parser.parsed("output_filename")) {
        auto output_filename = parser.get<std::string>("output_filename");
        std::cerr << "\tSaving data structure to disk...\n";
        essentials::save(locpres_mphf, output_filename.c_str());
        std::cerr << "\tDONE\n";
    }

    if (check) {
        std::cerr << "Checking\n";
        if (parser.parsed("output_filename")) {
            mphf_alt loaded;
            [[maybe_unused]] uint64_t num_bytes_read = essentials::load(loaded, parser.get<std::string>("output_filename").c_str());
            std::cerr << "[Info] Loaded " << num_bytes_read * 8 << " bits\n";
            pthash::bit_vector_builder population(loaded.get_kmer_count());  // bitvector for checking perfection and minimality
            if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
                std::cerr << "Unable to open the input file a second time" << input_filename << "\n";
                return 2;
            }
            seq = kseq_init(fp);
            while (check && kseq_read(seq) >= 0) 
            {
                std::string contig = std::string(seq->seq.s);
                check = check_collisions(loaded, contig, canonical, population);
                if (check) check = check_streaming_correctness(loaded, contig, canonical);
                // std::cerr << std::endl;
            }
            if (seq) kseq_destroy(seq);
            gzclose(fp);
            check = check && check_perfection(loaded, population);
        } else {
            pthash::bit_vector_builder population(locpres_mphf.get_kmer_count());
            if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
                std::cerr << "Unable to open the input file a second time" << input_filename << "\n";
                return 2;
            }
            seq = kseq_init(fp);
            while (check && kseq_read(seq) >= 0) 
            {
                std::string contig = std::string(seq->seq.s);
                check = check_collisions(locpres_mphf, contig, canonical, population);
                if (check) check = check_streaming_correctness(locpres_mphf, contig, canonical);
                // std::cerr << std::endl;
            }
            if (seq) kseq_destroy(seq);
            gzclose(fp);
            check = check && check_perfection(locpres_mphf, population);
        }
    }

    if (verbose) {
        std::cerr << "Statistics:\n";
        locpres_mphf.print_statistics();
    }

    if (total_contigs >= 1) --total_contigs;
    std::cout << input_filename << "," 
              << static_cast<uint32_t>(k) << "," 
              << static_cast<uint32_t>(m) << ","
              << static_cast<double>(total_colliding_minimizers) / total_distinct_minimizers << ","
              << 2.0/((k-m+1) + 1) << ","
              << static_cast<double>(total_minimizers) / total_kmers << ","
              << static_cast<double>(total_contigs) / total_kmers << ","
              << static_cast<double>(locpres_mphf.num_bits()) / locpres_mphf.get_kmer_count();
    std::cout << "\n";
    
    return 0;
}