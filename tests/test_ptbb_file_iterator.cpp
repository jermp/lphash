extern "C" {
#include "../external/kseq.h"
}
#include <zlib.h>
#include "../main/ptbb/ptbb.hpp"

using namespace lphash;

int main(int argc, char* argv[]) {
    if (argc < 3) return 1;
    std::string input_file(argv[1]);
    uint64_t k = atoll(argv[2]);
    ptbb::ptbb_file_itr kmer_itr(argv[1], k);
    ptbb::ptbb_file_itr kmer_end;
    std::size_t total_kmers = 0;
    for (; kmer_itr != kmer_end; ++kmer_itr) {
        ++total_kmers;
    }
    std::cerr << "total k-mers = " << total_kmers << std::endl;
    return 0;
}