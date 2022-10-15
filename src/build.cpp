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
    uint8_t k, m, nthreads;
    uint64_t mm_seed;
    std::string tmp_dirname;
    double c;
    bool canonical;
    bool in_memory;
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
    if (parser.parsed("in-memory")) in_memory = parser.get<bool>("in-memory");
    else in_memory = false;
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
    std::vector<mm_triplet_t> minimizers;
    std::cerr << "Part 1: file reading and info gathering\n";
    
    fp = NULL;
    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
        std::cerr << "Unable to open the input file " << input_filename << "\n";
        return 2;
    }
    seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        std::string contig = std::string(seq->seq.s);  // we lose a little bit of efficiency here
        auto n = minimizer::from_string<hash64>(contig, k, m, mm_seed, canonical, minimizers);  // non-canonical minimizers for now
        total_kmers += n;
        ++total_contigs;
        check_total_kmers += contig.length() - k + 1;
    }
    if (seq) kseq_destroy(seq);
    gzclose(fp);

    total_minimizers = minimizers.size();

    assert(total_kmers == check_total_kmers);

    std::cerr << "Part 2: build MPHF\n";
    mphf locpres_mphf(k, m, mm_seed, total_kmers, c, nthreads, in_memory, tmp_dirname, verbose);
    locpres_mphf.build_minimizers_mphf(minimizers);
    
    std::cerr << "Part 3: build fallback MPHF\n";
    // std::sort(colliding_minimizers.begin(), colliding_minimizers.end());
    { // garbage collector for unbucketable_kmers
        auto colliding_minimizers = locpres_mphf.build_inverted_index(minimizers);
        auto collect = []([[maybe_unused]] decltype(minimizers) to_collect) {};
        collect(std::move(minimizers));
        total_distinct_minimizers = locpres_mphf.get_minimizer_L0();
        total_colliding_minimizers = colliding_minimizers.size();
        std::vector<kmer_t> unbucketable_kmers;
        if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
            std::cerr << "Unable to open the input file a second time" << input_filename << "\n";
            return 2;
        }
        seq = kseq_init(fp);
        while (kseq_read(seq) >= 0) {
            std::string contig = std::string(seq->seq.s);
            minimizer::get_colliding_kmers<hash64>(contig, k, m, mm_seed, canonical, colliding_minimizers, unbucketable_kmers);
        }
        if (seq) kseq_destroy(seq);
        gzclose(fp);
        locpres_mphf.build_fallback_mphf(unbucketable_kmers);
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
            mphf loaded;
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