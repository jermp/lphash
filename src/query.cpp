extern "C" {
#include "../external/kseq.h"
}
#include <zlib.h>
#include "../include/partitioned_mphf.hpp"
#include "../include/unpartitioned_mphf.hpp"
#include "../include/query.hpp"
#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"

KSEQ_INIT(gzFile, gzread)

namespace lphash {

cmd_line_parser::parser get_query_parser(int argc, char* argv[]); // forward declaration

int query_partitioned_main(int argc, char* argv[]) {
    gzFile fp;
    kseq_t* seq;
    std::size_t total_kmers = 0;
    std::size_t total_dumb_kmers = 0;
    std::size_t total_aggregated_kmers = 0;
    std::size_t total_aggregated_dumb_kmers = 0;
    std::size_t total_time = 0;
    std::size_t total_dumb_time = 0;
    std::size_t total_aggregated_time = 0;
    std::size_t total_aggregated_dumb_time = 0;
    bool canonical = false;
    cmd_line_parser::parser parser(argc, argv);

    try {
        parser = get_query_parser(argc, argv);
    } catch (const ParseError& e) { return 1; }

    mphf hf;
    std::string mphf_filename = parser.get<std::string>("mphf");
    [[maybe_unused]] uint64_t num_bytes_read = essentials::load(hf, mphf_filename.c_str());
    std::string input_filename = parser.get<std::string>("input_filename");
    if (parser.parsed("canonical_parsing")) canonical = parser.get<bool>("canonical_parsing");

    essentials::timer<std::chrono::high_resolution_clock, std::chrono::microseconds> t;
    fp = NULL;
    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
        std::cerr << "Unable to open the input file " << input_filename << "\n";
        return 2;
    }
    seq = kseq_init(fp);
    t.start();
    while (kseq_read(seq) >= 0) {
        auto n = hf.barebone_streaming_query(seq->seq.s, seq->seq.l, canonical);
        total_kmers += n;
    }
    t.stop();
    if (seq) kseq_destroy(seq);
    gzclose(fp);
    total_time = t.elapsed();

    t.reset();

    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
        std::cerr << "Unable to open the input file " << input_filename << "\n";
        return 2;
    }
    seq = kseq_init(fp);
    t.start();
    while (kseq_read(seq) >= 0) {
        auto n = hf.barebone_dumb_query(seq->seq.s, seq->seq.l, canonical);
        total_dumb_kmers += n;
    }
    t.stop();
    if (seq) kseq_destroy(seq);
    gzclose(fp);
    total_dumb_time = t.elapsed();

    t.reset();

    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
        std::cerr << "Unable to open the input file " << input_filename << "\n";
        return 2;
    }
    seq = kseq_init(fp);
    t.start();
    while (kseq_read(seq) >= 0) {
        auto hashes = hf(seq->seq.s, seq->seq.l, canonical);
        total_aggregated_kmers += hashes.size();
        essentials::do_not_optimize_away(hashes.data());
    }
    t.stop();
    if (seq) kseq_destroy(seq);
    gzclose(fp);
    total_aggregated_time = t.elapsed();

    t.reset();

    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
        std::cerr << "Unable to open the input file " << input_filename << "\n";
        return 2;
    }
    seq = kseq_init(fp);
    t.start();
    while (kseq_read(seq) >= 0) {
        auto hashes = hf(seq->seq.s, seq->seq.l, canonical, false);
        total_aggregated_dumb_kmers += hashes.size();
        essentials::do_not_optimize_away(hashes.data());
    }
    t.stop();
    if (seq) kseq_destroy(seq);
    gzclose(fp);
    total_aggregated_dumb_time = t.elapsed();

    std::cout << input_filename << "," << mphf_filename << "," << total_kmers << ","
              << static_cast<double>(total_time * 1000) / total_kmers << ","
              << static_cast<double>(total_dumb_time * 1000) / total_dumb_kmers << ","
              << static_cast<double>(total_aggregated_time * 1000) / total_aggregated_kmers << ","
              << static_cast<double>(total_aggregated_dumb_time * 1000) /
                     total_aggregated_dumb_kmers
              << std::endl;

    // std::cerr << "number of k-mers:"
    //           << " " << total_kmers << " " << total_dumb_kmers << " " << total_aggregated_kmers
    //           << " " << total_aggregated_dumb_kmers << std::endl;
    return 0;
}

int query_unpartitioned_main(int argc, char* argv[]) {
    gzFile fp;
    kseq_t* seq;
    std::size_t total_kmers = 0;
    std::size_t total_dumb_kmers = 0;
    std::size_t total_aggregated_kmers = 0;
    std::size_t total_aggregated_dumb_kmers = 0;
    std::size_t total_time = 0;
    std::size_t total_dumb_time = 0;
    std::size_t total_aggregated_time = 0;
    std::size_t total_aggregated_dumb_time = 0;
    bool canonical = false;
    cmd_line_parser::parser parser(argc, argv);

    try {
        parser = get_query_parser(argc, argv);
    } catch (const ParseError& e) { return 1; }

    mphf_alt hf;
    std::string mphf_filename = parser.get<std::string>("mphf");
    [[maybe_unused]] uint64_t num_bytes_read = essentials::load(hf, mphf_filename.c_str());
    std::string input_filename = parser.get<std::string>("input_filename");
    if (parser.parsed("canonical_parsing")) canonical = parser.get<bool>("canonical_parsing");

    essentials::timer<std::chrono::high_resolution_clock, std::chrono::microseconds> t;
    fp = NULL;
    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
        std::cerr << "Unable to open the input file " << input_filename << "\n";
        return 2;
    }
    seq = kseq_init(fp);
    t.start();
    while (kseq_read(seq) >= 0) {
        auto n = hf.barebone_streaming_query(seq->seq.s, seq->seq.l, canonical);
        total_kmers += n;
    }
    t.stop();
    if (seq) kseq_destroy(seq);
    gzclose(fp);
    total_time = t.elapsed();

    t.reset();

    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
        std::cerr << "Unable to open the input file " << input_filename << "\n";
        return 2;
    }
    seq = kseq_init(fp);
    t.start();
    while (kseq_read(seq) >= 0) {
        auto n = hf.barebone_dumb_query(seq->seq.s, seq->seq.l, canonical);
        total_dumb_kmers += n;
    }
    t.stop();
    if (seq) kseq_destroy(seq);
    gzclose(fp);
    total_dumb_time = t.elapsed();

    t.reset();

    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
        std::cerr << "Unable to open the input file " << input_filename << "\n";
        return 2;
    }
    seq = kseq_init(fp);
    t.start();
    while (kseq_read(seq) >= 0) {
        auto hashes = hf(seq->seq.s, seq->seq.l, canonical);
        total_aggregated_kmers += hashes.size();
        essentials::do_not_optimize_away(hashes.data());
    }
    t.stop();
    if (seq) kseq_destroy(seq);
    gzclose(fp);
    total_aggregated_time = t.elapsed();

    t.reset();

    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
        std::cerr << "Unable to open the input file " << input_filename << "\n";
        return 2;
    }
    seq = kseq_init(fp);
    t.start();
    while (kseq_read(seq) >= 0) {
        auto hashes = hf(seq->seq.s, seq->seq.l, canonical, false);
        total_aggregated_dumb_kmers += hashes.size();
        essentials::do_not_optimize_away(hashes.data());
    }
    t.stop();
    if (seq) kseq_destroy(seq);
    gzclose(fp);
    total_aggregated_dumb_time = t.elapsed();

    std::cout << input_filename << "," << mphf_filename << "," << total_kmers << ","
              << static_cast<double>(total_time * 1000) / total_kmers << ","
              << static_cast<double>(total_dumb_time * 1000) / total_dumb_kmers << ","
              << static_cast<double>(total_aggregated_time * 1000) / total_aggregated_kmers << ","
              << static_cast<double>(total_aggregated_dumb_time * 1000) /
                     total_aggregated_dumb_kmers
              << std::endl;
    return 0;
}

cmd_line_parser::parser get_query_parser(int argc, char* argv[]) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("mphf", "lphash minimal perfect hash function saved on disk\n");
    parser.add("input_filename",
               "Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not:\n");
    if (!parser.parse()) throw ParseError();
    return parser;
}

}