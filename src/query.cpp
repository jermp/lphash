extern "C" {
#include "../external/kseq.h"
}
#include <zlib.h>

#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../external/pthash/external/essentials/include/essentials.hpp"
#include "../include/util.hpp"

KSEQ_INIT(gzFile, gzread)

namespace lphash {

cmd_line_parser::parser get_query_parser(int argc, char* argv[]) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("mphf", "LP-Hash MPHF saved on disk.\n", "-i", true);
    parser.add("query_filename",
               "Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not.",
               "-q", true);
    if (!parser.parse()) throw ParseError();
    return parser;
}

template <typename MPHF>
int query(int argc, char* argv[]) {
    gzFile fp;
    kseq_t* seq;

    cmd_line_parser::parser parser(argc, argv);

    try {
        parser = get_query_parser(argc, argv);
    } catch (const ParseError& e) { return 1; }

    MPHF hf;
    std::string mphf_filename = parser.get<std::string>("mphf");
    essentials::load(hf, mphf_filename.c_str());
    std::string query_filename = parser.get<std::string>("query_filename");

    essentials::timer<std::chrono::high_resolution_clock, std::chrono::microseconds> t;
    fp = NULL;

    if ((fp = gzopen(query_filename.c_str(), "r")) == NULL) {
        std::cerr << "Unable to open the input file " << query_filename << "\n";
        return 2;
    }
    seq = kseq_init(fp);
    t.start();
    constexpr bool streaming_on = true;
    uint64_t total_kmers_streaming_on = 0;
    while (kseq_read(seq) >= 0) {
        auto hashes = hf(seq->seq.s, seq->seq.l, streaming_on);
        total_kmers_streaming_on += hashes.size();
        essentials::do_not_optimize_away(hashes.data());
    }
    t.stop();
    if (seq) kseq_destroy(seq);
    gzclose(fp);
    uint64_t time_streaming_on = t.elapsed();

    t.reset();

    if ((fp = gzopen(query_filename.c_str(), "r")) == NULL) {
        std::cerr << "Unable to open the input file " << query_filename << "\n";
        return 2;
    }
    seq = kseq_init(fp);
    t.start();
    constexpr bool streaming_off = false;
    uint64_t total_kmers_streaming_off = 0;
    while (kseq_read(seq) >= 0) {
        auto hashes = hf(seq->seq.s, seq->seq.l, streaming_off);
        total_kmers_streaming_off += hashes.size();
        essentials::do_not_optimize_away(hashes.data());
    }
    t.stop();
    if (seq) kseq_destroy(seq);
    gzclose(fp);
    uint64_t time_streaming_off = t.elapsed();

    assert(total_kmers_streaming_on == total_kmers_streaming_off);

    std::cout << query_filename << "," << mphf_filename << "," << total_kmers_streaming_on << ","
              << static_cast<double>(time_streaming_on * 1000) / total_kmers_streaming_on << ","
              << static_cast<double>(time_streaming_off * 1000) / total_kmers_streaming_off
              << std::endl;
    return 0;
}

}  // namespace lphash