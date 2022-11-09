extern "C" {
#include "../external/kseq.h"
}

#include "../include/constants.hpp"
#include "../include/ptbb_file_itr.hpp"
#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../external/BooPHF.hpp"

KSEQ_INIT(gzFile, gzread)

using namespace lphash;

int main(int argc, char* argv[]) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("input_filename",
               "Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not. "
               "All k-mers must have been seen by the MPHFs.");
    parser.add("k", "K-mer length (must be <= " + std::to_string(constants::max_k) + ").");
    parser.add("pthash_filename", "Output file name where the pthash mphf will be serialized.",
               "-p", false);
    parser.add("bbhash_filename", "Output file name where the BBHash mphf will be serialized.",
               "-b", false);
    if (!parser.parse()) return 1;

    std::string input_filename = parser.get<std::string>("input_filename");
    uint64_t k = parser.get<uint64_t>("k");

    std::cout << input_filename << "," << k;

    pthash_mphf_t kmer_order;
    if (parser.parsed("pthash_filename")) {
        std::string pthash_filename = parser.get<std::string>("pthash_filename");
        essentials::load(kmer_order, pthash_filename.c_str());
    }
    auto bphf = boomphf::mphf<kmer_t, other::BBHasher<kmer_t>>();
    if (parser.parsed("bbhash_filename")) {
        std::string bbhash_filename = parser.get<std::string>("bbhash_filename");
        std::ifstream bbh_strm(bbhash_filename, std::ios::binary);
        bphf.load(bbh_strm);
    }
    /*
    std::size_t total_kmers = 0;
    std::size_t pt_ns = 0, bb_ns = 0;
    {
    other::ptbb_file_itr kmer_itr(input_filename, k);
    // kseq_t* itr_guts = reinterpret_cast<kseq_t*>(kmer_itr.memory_management());
    other::ptbb_file_itr kmer_end;
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::nanoseconds> pt_timer;
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::nanoseconds> bb_timer;
    for (; kmer_itr != kmer_end; ++kmer_itr) {
            kmer_t kmer = *kmer_itr;
            if (parser.parsed("pthash_filename")) {
                    pt_timer.start();
                    [[maybe_unused]] auto hval = kmer_order(kmer);
                    pt_timer.stop();
                    pt_ns += pt_timer.elapsed();
                    pt_timer.reset();
                    essentials::do_not_optimize_away(hval);
            }
            if (parser.parsed("bbhash_filename")) {
                    bb_timer.start();
                    [[maybe_unused]] auto hval = bphf.lookup(kmer);
                    bb_timer.stop();
                    bb_ns += bb_timer.elapsed();
                    bb_timer.reset();
                    essentials::do_not_optimize_away(hval);
            }
            ++total_kmers;
    }
    }// kseq_destroy(itr_guts);
    */
    std::size_t pt_us = 0, bb_us = 0;
    std::size_t pt_total_kmers = 1, bb_total_kmers = 1;
    if (parser.parsed("pthash_filename")) {
        pt_total_kmers = 0;
        other::ptbb_file_itr kmer_itr(input_filename, k);
        other::ptbb_file_itr kmer_end;
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::microseconds> pt_timer;
        pt_timer.start();
        for (; kmer_itr != kmer_end; ++kmer_itr) {
            [[maybe_unused]] auto hval = kmer_order(*kmer_itr);
            essentials::do_not_optimize_away(hval);
            ++pt_total_kmers;
        }
        pt_timer.stop();
        pt_us = pt_timer.elapsed();
    }
    if (parser.parsed("bbhash_filename")) {
        bb_total_kmers = 0;
        other::ptbb_file_itr kmer_itr(input_filename, k);
        other::ptbb_file_itr kmer_end;
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::microseconds> bb_timer;
        bb_timer.start();
        for (; kmer_itr != kmer_end; ++kmer_itr) {
            [[maybe_unused]] auto hval = bphf.lookup(*kmer_itr);
            essentials::do_not_optimize_away(hval);
            ++bb_total_kmers;
        }
        bb_timer.stop();
        bb_us = bb_timer.elapsed();
    }

    if (parser.parsed("pthash_filename")) {
        std::cout << "," << parser.get<std::string>("pthash_filename") << ","
                  << static_cast<double>(pt_us * 1000) / pt_total_kmers;
    } else {
        std::cout << ",,";
    }
    if (parser.parsed("bbhash_filename")) {
        std::cout << "," << parser.get<std::string>("bbhash_filename") << ","
                  << static_cast<double>(bb_us * 1000) / bb_total_kmers;
    } else {
        std::cout << ",,";
    }
    if (parser.parsed("pthash_filename") && parser.parsed("bbhash_filename"))
        assert(pt_total_kmers == bb_total_kmers);
    std::cout << "\n";
}