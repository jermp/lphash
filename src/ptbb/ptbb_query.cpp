#include <iostream>

#include "../../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "ptbb.hpp"

using namespace lphash;

int main(int argc, char* argv[]) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("input_filename",
               "Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not.");
    parser.add("k", "K-mer length (must be <= " + std::to_string(constants::max_k) + ").");
    parser.add("pthash_filename", "PTHash MPHF",
               "-p", false);
    parser.add("bbhash_filename", "BBHash MPHF",
               "-b", false);
    // parser.add("sichash_filename", "SicHash MPHF", "-s", false); Does not work since SicHash doesn't support saving
    if (!parser.parse()) return 1;

    std::string input_filename = parser.get<std::string>("input_filename");
    uint64_t k = parser.get<uint64_t>("k");

    std::cout << input_filename << "," << k;

    ptbb::pthash_mphf_t pthash_mphf;
    if (parser.parsed("pthash_filename")) {
        std::string pthash_filename = parser.get<std::string>("pthash_filename");
        essentials::load(pthash_mphf, pthash_filename.c_str());
    }

    ptbb::bbhash_mphf_t bbhash_mphf;
    if (parser.parsed("bbhash_filename")) {
        std::string bbhash_filename = parser.get<std::string>("bbhash_filename");
        std::ifstream bbh_strm(bbhash_filename, std::ios::binary);
        bbhash_mphf.load(bbh_strm);
    }

    uint64_t pt_us = 0, bb_us = 0;
    uint64_t pt_total_kmers = 1, bb_total_kmers = 1;
    if (parser.parsed("pthash_filename")) {
        pt_total_kmers = 0;
        ptbb::ptbb_file_itr kmer_itr(input_filename, k);
        ptbb::ptbb_file_itr kmer_end;
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::microseconds> pt_timer;
        pt_timer.start();
        for (; kmer_itr != kmer_end; ++kmer_itr) {
            auto hval = pthash_mphf(*kmer_itr);
            essentials::do_not_optimize_away(hval);
            ++pt_total_kmers;
        }
        pt_timer.stop();
        pt_us = pt_timer.elapsed();
    }
    if (parser.parsed("bbhash_filename")) {
        bb_total_kmers = 0;
        ptbb::ptbb_file_itr kmer_itr(input_filename, k);
        ptbb::ptbb_file_itr kmer_end;
        essentials::timer<std::chrono::high_resolution_clock, std::chrono::microseconds> bb_timer;
        bb_timer.start();
        for (; kmer_itr != kmer_end; ++kmer_itr) {
            auto hval = bbhash_mphf.lookup(*kmer_itr);
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
    std::cout << "\n";
}