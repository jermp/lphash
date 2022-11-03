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
    } catch (const ParseError& e) {
        return 1;
    }

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
}

/*
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
    gzclose(fp);
*/