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
    cmd_line_parser::parser parser = get_query_parser(argc, argv);
    mphf hf;
    std::string mphf_filename = parser.get<std::string>("mphf");
    std::cerr << "Loading mphf: " << mphf_filename << " ...\n";
    [[maybe_unused]] uint64_t num_bytes_read = essentials::load(hf, mphf_filename.c_str());
    std::cerr << "Loading DONE\n";
    std::cerr << hf << "\n";
    std::string input_filename = parser.get<std::string>("input_filename");

    std::cerr << "Input file = " << input_filename << std::endl;

    fp = NULL;
    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
        std::cerr << "Unable to open the input file " << input_filename << "\n";
        return 2;
    }
    
    seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        std::string contig = std::string(seq->seq.s);  // we lose a little bit of efficiency here
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::microseconds> t;
        t.start();
        auto fast_hashes = hf(contig);
        t.stop();
        total_time += t.elapsed();
        total_kmers += fast_hashes.size();
    }
    if (seq) kseq_destroy(seq);
    std::cout << (total_time * 1000) / total_kmers << " ns/kmer" << std::endl;
}