#include <iostream>

#include "../include/constants.hpp"
#include "ptbb.hpp"

#include <argparse/argparse.hpp>

using namespace lphash;

int main(int argc, char* argv[]) {
    argparse::ArgumentParser parser(argv[0]);
    parser.add_argument("-q", "--query-filename")
        .help("Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not.")
        .required();
    parser.add_argument("-k")
        .help("K-mer length (must be <= " + std::to_string(constants::max_k) + ").")
        .scan<'u', uint64_t>()
        .required();
    parser.add_argument("-p", "--pthash-filename")
        .help("PTHash MPHF");
    parser.add_argument("-b", "--bbhash-filename")
        .help("BBHash MPHF");
    
    try {
        parser.parse_args(argc, argv);
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    std::string input_filename = parser.get<std::string>("-q");
    uint64_t k = parser.get<uint64_t>("-k");

    std::cout << input_filename << "," << k;

    ptbb::pthash_mphf_t pthash_mphf;
    if (parser.is_used("--pthash-filename")) {
        std::string pthash_filename = parser.get<std::string>("--pthash-filename");
        essentials::load(pthash_mphf, pthash_filename.c_str());
    }

    ptbb::bbhash_mphf_t bbhash_mphf;
    if (parser.is_used("--bbhash-filename")) {
        std::string bbhash_filename = parser.get<std::string>("--bbhash-filename");
        std::ifstream bbh_strm(bbhash_filename, std::ios::binary);
        bbhash_mphf.load(bbh_strm);
    }

    uint64_t pt_us = 0, bb_us = 0;
    uint64_t pt_total_kmers = 1, bb_total_kmers = 1;
    if (parser.is_used("--pthash-filename")) {
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
    if (parser.is_used("--bbhash-filename")) {
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

    if (parser.is_used("--pthash-filename")) {
        std::cout << "," << parser.get<std::string>("--pthash-filename") << ","
                  << static_cast<double>(pt_us * 1000) / pt_total_kmers;
    } else {
        std::cout << ",,";
    }
    if (parser.is_used("--bbhash-filename")) {
        std::cout << "," << parser.get<std::string>("--bbhash-filename") << ","
                  << static_cast<double>(bb_us * 1000) / bb_total_kmers;
    } else {
        std::cout << ",,";
    }
    std::cout << "\n";
}