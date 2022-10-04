extern "C" {
#include "../include/kseq.h"
}
#include <zlib.h>
#include "../include/mphf.hpp"
#include "../include/parser_build.hpp"

KSEQ_INIT(gzFile, gzread)

using namespace lphash;

int main(int argc, char* argv[]) {
    gzFile fp;
    kseq_t* seq;
    std::size_t total_kmers = 0;
    std::size_t total_time = 0;
    std::size_t total_dumb_time = 0;
    bool canonical = false;
    cmd_line_parser::parser parser = get_query_parser(argc, argv);
    mphf hf;
    std::string mphf_filename = parser.get<std::string>("mphf");
    // std::cerr << "Loading mphf: " << mphf_filename << " ...\n";
    [[maybe_unused]] uint64_t num_bytes_read = essentials::load(hf, mphf_filename.c_str());
    // std::cerr << "Loading DONE\n";
    // std::cerr << hf << "\n";
    std::string input_filename = parser.get<std::string>("input_filename");

    // std::cerr << "Input file = " << input_filename << std::endl;

    fp = NULL;
    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
        std::cerr << "Unable to open the input file " << input_filename << "\n";
        return 2;
    }
    
    if (parser.parsed("canonical_parsing")) canonical = parser.get<bool>("canonical_parsing");

    seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        std::string contig = std::string(seq->seq.s);  // we lose a little bit of efficiency here
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::nanoseconds> t;
        t.reset();
        t.start();
        auto n = hf.barebone_streaming_query(contig, canonical);
        t.stop();
        total_time += t.elapsed();
        total_kmers += n;
        t.reset();
        t.start();
        auto dumb_hashes = hf.dumb_evaluate(contig, false);
        t.stop();
        total_dumb_time += t.elapsed();
    }
    if (seq) kseq_destroy(seq);
    std::cout << input_filename << "," << mphf_filename << "," << static_cast<double>(total_time) / total_kmers << "," << static_cast<double>(total_dumb_time) / total_kmers << std::endl;
}