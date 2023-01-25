#include "../include/constants.hpp"
#include "../include/parser_build.hpp"

#include <chrono>

namespace lphash {

typedef std::chrono::high_resolution_clock clock_type;

template <typename MPHF>
int build(int argc, char* argv[]) {
    configuration config;
    try {
        parse_build_config(argc, argv, config);
    } catch (const ParseError& e) {
        std::cerr << e.what() << std::endl;
        return 2;
    } catch (const OptionError& e) {
        std::cerr << e.what() << std::endl;
        return 3;
    }

    auto start = clock_type::now();
    MPHF f;
    f.build(config, std::cout);
    if (config.output_filename != "") {
        if (config.verbose) std::cerr << "Saving data structure to disk...  ";
        essentials::save(f, config.output_filename.c_str());
        if (config.verbose) std::cerr << "DONE\n";
    }
    auto stop = clock_type::now();
    std::cout << "function built in "
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

template <typename MPHF>
void check(MPHF const& f, configuration const& config) {
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

}  // namespace lphash