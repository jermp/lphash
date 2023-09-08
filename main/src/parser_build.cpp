#include "../../lib/include/constants.hpp"
#include "../include/parser_build.hpp"
#include "../../external/pthash/external/cmd_line_parser/include/parser.hpp"

namespace lphash {

cmd_line_parser::parser get_build_parser(int argc, char* argv[]) {
    cmd_line_parser::parser parser(argc, argv);

    /* mandatory arguments */
    parser.add("input_filename",
               "Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not:\n"
               "\t- without duplicate nor invalid kmers\n"
               "\t- one DNA sequence per line.\n"
               "\tFor example, it could be the de Bruijn graph topology output by BCALM.",
               "-i", true);
    parser.add("k", "K-mer length (must be <= " + std::to_string(constants::max_k) + ").", "-k",
               true);
    parser.add("m", "Minimizer length (must be < k and <= 32).", "-m", true);

    /* optional arguments */
    parser.add("seed",
               "Seed for minimizer computation (default is " +
                   std::to_string(constants::default_seed) + ").",
               "-s", false);
    parser.add("threads",
               "Number of threads for pthash (default is " +
                   std::to_string(constants::default_num_threads) + ").",
               "-t", false);
    parser.add("output_filename",
               "Output file name where the data structure will be serialized (no files generated "
               "by default).",
               "-o", false);
    parser.add(
        "tmp_dirname",
        "Temporary directory used for construction in external memory (default is directory '" +
            constants::default_tmp_dirname + "').",
        "-d", false);
    parser.add(
        "c",
        "A (floating point) constant that trades construction speed for space effectiveness of "
        "minimal perfect hashing. \n"
        "\tA reasonable value lies between 3.0 and 10.0 (default is " +
            std::to_string(constants::c).substr(0, std::to_string(constants::c).find(".") + 2 + 1) +
            ").",
        "-c", false);
    parser.add("max-memory",
               "Maximum internal memory [GB] for building (8GB by default). Use external memory if "
               "needed.",
               "--max-memory", false);
    parser.add("check", "Check correctness after construction (disabled by default).", "--check",
               false, true);
    parser.add("verbose", "Verbose output during construction (disabled by default).", "--verbose",
               false, true);

    if (!parser.parse()) throw ParseError();
    return parser;
}

void parse_build_config(int argc, char* argv[], configuration& config) {
    auto parser = get_build_parser(argc, argv);

    config.input_filename = parser.get<std::string>("input_filename");
    config.k = parser.get<uint64_t>("k");
    if (config.k > constants::max_k)
        throw OptionError("k cannot be larger than " + std::to_string(constants::max_k));
    config.m = parser.get<uint64_t>("m");
    if (config.m > config.k) throw OptionError("m cannot be larger than k");
    if (parser.parsed("output_filename"))
        config.output_filename = parser.get<std::string>("output_filename");

    if (parser.parsed("seed")) config.mm_seed = parser.get<uint64_t>("seed");
    if (parser.parsed("threads")) config.num_threads = parser.get<uint64_t>("threads");
    if (parser.parsed("tmp_dirname")) {
        config.tmp_dirname = parser.get<std::string>("tmp_dirname");
        essentials::create_directory(config.tmp_dirname);
    }

    if (parser.parsed("c")) {
        config.c = parser.get<double>("c");
        if (config.c > 10.0 || config.c < 3.0) throw OptionError("3.0 <= c <= 10.0");
    }

    if (parser.parsed("max-memory")) {
        config.max_memory = parser.get<uint64_t>("max-memory");
        if (config.max_memory > 255)
            throw OptionError("The maximum allowed amount of ram is 255GB");
    }

    if (parser.parsed("check")) config.check = parser.get<bool>("check");
    if (parser.parsed("verbose")) config.verbose = parser.get<bool>("verbose");
}

}  // namespace lphash
