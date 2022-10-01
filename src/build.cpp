extern "C" {
#include "../include/kseq.h"
}
#include <zlib.h>
#include "../include/constants.hpp"
#include "../include/parser_build.hpp"
#include "minimizer.hpp"
#include "../include/mphf.hpp"

#include "../include/quartet_wtree.hpp"
#include "../include/prettyprint.hpp"

using namespace lphash;

KSEQ_INIT(gzFile, gzread)

int main(int argc, char* argv[]) {
    gzFile fp;
    kseq_t* seq;
    uint8_t k, m;
    uint64_t mm_seed;
    std::string tmp_dirname;
    bool verbose;
    bool check;

    std::size_t total_kmers;

    cmd_line_parser::parser parser = get_build_parser(argc, argv);
    std::string input_filename = parser.get<std::string>("input_filename");
    k = static_cast<uint8_t>(parser.get<uint32_t>("k"));
    m = static_cast<uint8_t>(parser.get<uint32_t>("m"));
    if (parser.parsed("seed")) mm_seed = parser.get<uint64_t>("seed");
    else mm_seed = 42;
    if (parser.parsed("tmp_dirname")) tmp_dirname = parser.get<std::string>("tmp_dirname");
    else tmp_dirname = "";
    if (parser.parsed("check")) check = parser.get<bool>("check");
    else check = false;
    if (parser.parsed("verbose")) verbose = parser.get<bool>("verbose");
    else verbose = false;

    std::cerr << uint32_t(k) << " " << uint32_t(m);
    std::cerr << "mm_seed = " << mm_seed;
    std::cerr << "\n";

    total_kmers = 0;
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
        auto n = minimizer::from_string<hash64>(contig, k, m, mm_seed, false, minimizers);  // not canonical minimizers for now
        total_kmers += n;
    }
    if (seq) kseq_destroy(seq);

    std::cerr << "Part 2: build MPHF\n";
    mphf locpres_mphf(k, m, mm_seed, total_kmers, 1, tmp_dirname, verbose);
    locpres_mphf.build_minimizers_mphf(minimizers);
    auto colliding_minimizers = locpres_mphf.build_inverted_index(minimizers);
    
    std::cerr << "Part 3: build fallback MPHF\n";
    std::sort(colliding_minimizers.begin(), colliding_minimizers.end());  // FIXME sort minimizers for fast search -> find better alternative (hash table)
    { // garbage collector for unbucketable_kmers
        std::vector<kmer_t> unbucketable_kmers;
        unbucketable_kmers.reserve(total_kmers);  // worst case scenario
        if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {  // reopen input file in order to find ambiguous minimizers and their k-mers
            std::cerr << "Unable to open the input file a second time" << input_filename << "\n";
            return 2;
        }
        seq = kseq_init(fp);
        while (kseq_read(seq) >= 0) {
            std::string contig = std::string(seq->seq.s);  // we lose a little bit of efficiency here
            minimizer::get_colliding_kmers<hash64, kmer_t>(contig, k, m, mm_seed, false, colliding_minimizers, unbucketable_kmers);
        }
        if (seq) kseq_destroy(seq);
        locpres_mphf.build_fallback_mphf(unbucketable_kmers);
    }

    if (parser.parsed("output_filename")) {
        auto output_filename = parser.get<std::string>("output_filename");
        std::cerr << "saving data structure to disk\n";
        essentials::save(locpres_mphf, output_filename.c_str());
        std::cerr << "DONE\n";
    }

    if (check) {
        std::cerr << "Part 4: check\n";
        if (parser.parsed("output_filename")) {
            mphf loaded;
            uint64_t num_bytes_read = essentials::load(loaded, parser.get<std::string>("output_filename").c_str());
            locpres_mphf = loaded;
        }
        pthash::bit_vector_builder population(locpres_mphf.get_kmer_count());  // bitvector for checking perfection and minimality
        if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {  // reopen input stream once again
            std::cerr << "Unable to open the input file a second time" << input_filename << "\n";
            return 2;
        }
        seq = kseq_init(fp);
        while (check && kseq_read(seq) >= 0) 
        {
            std::string contig = std::string(seq->seq.s);  // we lose a little bit of efficiency here
            check = check_collisions(locpres_mphf, contig, population);
        }
        if (seq) kseq_destroy(seq);
        check = check_perfection(locpres_mphf, population);
    }
    std::cerr << "Statistics:\n";
    locpres_mphf.print_statistics();

    
    
    return 0;
}