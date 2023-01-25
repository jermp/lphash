#include "../include/constants.hpp"
#include "../include/partitioned_mphf.hpp"
#include "../include/unpartitioned_mphf.hpp"
#include "../include/parser_build.hpp"

#include <chrono>

namespace lphash {

typedef std::chrono::high_resolution_clock clock_type;

int build_partitioned_main(int argc, char* argv[]) {
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
    mphf f;
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

int build_unpartitioned_main(int argc, char* argv[]) {
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
    mphf_alt f;
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
        check_alt(f, config);
    }

    if (config.verbose) {
        std::cerr << "Statistics:\n";
        f.print_statistics();
    }

    return 0;
}

}  // namespace lphash