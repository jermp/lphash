#pragma once

#include "../../external/pthash/external/essentials/include/essentials.hpp"
#include "../../external/pthash/include/encoders/bit_vector.hpp"
#include "../include/constants.hpp"

#include <chrono>
#include <argparse/argparse.hpp>

namespace lphash {

argparse::ArgumentParser get_parser_build();

template <class MPHF>
typename MPHF::configuration get_configuration(const argparse::ArgumentParser& parser) {
    typename MPHF::configuration config;
    config.input_filename = parser.get<std::string>("-i");
    config.output_filename = parser.get<std::string>("-o");
    config.k = parser.get<uint64_t>("-k");
    config.m = parser.get<uint64_t>("-m");
    config.mm_seed = parser.get<uint64_t>("-s");
    config.pt_seed = parser.get<uint64_t>("--pthash_seed");
    config.c = parser.get<double>("-c");
    config.num_threads = parser.get<uint64_t>("-t");
    config.max_memory = parser.get<uint64_t>("--max-memory");
    config.tmp_dirname = parser.get<std::string>("--tmp-dir");
    config.check = parser.get<bool>("--check");
    config.verbose = parser.get<bool>("--verbose");

    if (config.k > constants::max_k) throw OptionError("k cannot be larger than " + std::to_string(constants::max_k));
    if (config.m > config.k) throw OptionError("m cannot be larger than k");
    if (config.tmp_dirname != ".") essentials::create_directory(config.tmp_dirname);
    if (config.c > 10.0 || config.c < 3.0) throw OptionError("3.0 <= c <= 10.0");
    // if (config.max_memory > 255) throw OptionError("The maximum allowed amount of ram is 255 GB");

    return config;
}

template <typename MPHF>
void check(MPHF const& f, typename MPHF::configuration const& config) {
    gzFile fp = nullptr;
    kseq_t* seq = nullptr;
    pthash::bit_vector_builder population(f.get_kmer_count());
    if ((fp = gzopen(config.input_filename.c_str(), "r")) == NULL)
        throw std::runtime_error("Unable to open input file " + config.input_filename +
                                 " for checking\n");
    seq = kseq_init(fp);
    bool good = true;
    while (good && kseq_read(seq) >= 0) {
        good = check_collisions(f, seq->seq.s, seq->seq.l, population);
        if (good) good = check_streaming_correctness(f, seq->seq.s, seq->seq.l);
    }
    if (seq) kseq_destroy(seq);
    gzclose(fp);
    if (good) check_perfection(f, population);
}

template <typename MPHF>
int build_main(const argparse::ArgumentParser& parser) {
    typedef std::chrono::high_resolution_clock clock_type;
    auto config = get_configuration<MPHF>(parser);
    // try {
    //     parse_build_config(argc, argv, config);
    // } catch (const ParseError& e) {
    //     std::cerr << e.what() << std::endl;
    //     return 2;
    // } catch (const OptionError& e) {
    //     std::cerr << e.what() << std::endl;
    //     return 3;
    // }

    auto start = clock_type::now();
    MPHF f;
    f.build(config, std::cout);
    if (config.output_filename != "") {
        if (config.verbose) std::cerr << "Saving data structure to disk...  ";
        essentials::save(f, config.output_filename.c_str());
        if (config.verbose) std::cerr << "DONE\n";
    }
    auto stop = clock_type::now();
    std::cerr << "function built in "
              << std::chrono::duration_cast<std::chrono::seconds>(stop - start).count() << " [sec]"
              << std::endl;

    if (config.check) {
        std::cerr << "Checking...\n";
        if (config.output_filename != "") {
            uint64_t num_bytes_read = essentials::load(f, config.output_filename.c_str());
            std::cerr << "[Info] Loaded " << num_bytes_read * 8 << " bits\n";
        }
        check(f, config);
    }

    if (config.verbose) {
        std::cerr << "Statistics:\n";
        f.print_statistics();
    }

    return 0;
}

} // namespace lphash