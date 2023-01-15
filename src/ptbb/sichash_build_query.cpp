#include <iostream>

#include "../../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "ptbb.hpp"

using namespace lphash;

/*
 * SicHash tests are implemented in this standalone file 
 * since ptbb works by saving/loading the hash table 
 * which are not supported by SicHash
 */

int main(int argc, char* argv[]) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("input_filename",
               "FASTA file used to build SicHash MPHF (.fa/fasta extension) compressed with gzip (.gz) or not:\n"
               "\t- without duplicate nor invalid kmers\n"
               "\t- one DNA sequence per line.\n"
               "\tFor example, it could be the de Bruijn graph topology output by BCALM.");
    parser.add("query_filename", "Query FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not");
    parser.add("k", "K-mer length (must be <= " + std::to_string(constants::max_k) + ").");
    parser.add("p1", "SicHash p1", "-p", false);
    parser.add("p2", "SicHash p2", "-P", false);
    // parser.add(
    //     "tmp_dirname",
    //     "Temporary directory used for construction in external memory. Default is directory '" +
    //         constants::default_tmp_dirname + "'.",
    //     "-d", false);
    // parser.add("threads", "Number of threads for pthash (default is 1).", "-t", false);
    // parser.add("verbose", "Verbose output during construction.", "--verbose", true);
    parser.add("check", "Check output", "--check", true);
    if (!parser.parse()) return 1;

    auto input_filename = parser.get<std::string>("input_filename");
    auto k = parser.get<uint64_t>("k");
    bool check = parser.get<bool>("check");
    std::string input_filename = parser.get<std::string>("input_filename");
    std::string query_filename = parser.get<std::string>("query_filename");

    sichash::SicHashConfig sic_config;
    double p1 = parser.get<double>("p1");
    double p2 = parser.get<double>("p2");
    sic_config.percentages(p1, p2);
    sic_config.silent = true;
    auto kmers = ptbb::load_kmers(input_filename, k);
    sichash::SicHash<true> sichash_mphf(kmers, sic_config);
    std::cout << input_filename << ","  << k << "," << sichash_mphf.spaceUsage() << "," << static_cast<double>(sichash_mphf.spaceUsage()) / kmers.size();

    if (check) {
        std::cerr << "Checking SicHash...";
        pthash::bit_vector_builder population(kmers.size());
        uint64_t check_total_kmers = 0;
        {
            for (auto kmer_itr = kmers.begin(); kmer_itr != kmers.end(); ++kmer_itr) {
                auto idx = sichash_mphf(*kmer_itr);
                if (idx >= kmers.size()) {
                    std::cerr << "[Error] out of bounds: " << idx << std::endl;
                    return 2;
                } else if (population.get(idx)) {
                    std::cerr << "[Error] collision" << std::endl;
                    return 2;
                } else
                    population.set(idx);
                ++check_total_kmers;
            }
        }
        assert(kmers.size() == check_total_kmers);
        for (uint64_t i = 0; i < kmers.size(); ++i) {
            if (!population.get(i)) {
                std::cerr << "[Error] hash is not perfect" << std::endl;
                return 2;
            }
        }
        std::cerr << "EVERYTHING OK\n";
    }

    kmers = ptbb::load_kmers(query_filename, k);
    essentials::timer<std::chrono::high_resolution_clock, std::chrono::microseconds> pt_timer;
    pt_timer.start();
    for (auto& kmer : kmers) {
        auto hval = sichash_mphf(kmer);
        essentials::do_not_optimize_away(hval);
    }
    pt_timer.stop();

    auto pt_us = pt_timer.elapsed();
    std::cout << "," << query_filename << "," << static_cast<double>(pt_us * 1000) / kmers.size() << "\n";
}