#pragma once

#include <exception>
#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"

namespace lphash {

class ParseError : public std::exception {
public:
    const char* what() const noexcept { return "Unable to parse the arguments\n"; }
};

cmd_line_parser::parser get_build_parser(int argc, char* argv[]) {
    cmd_line_parser::parser parser(argc, argv);

    /* mandatory arguments */
    parser.add("input_filename",
               "Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not:\n"
               "\t- without duplicate nor invalid kmers\n"
               "\t- one DNA sequence per line.\n"
               "\tFor example, it could be the de Bruijn graph topology output by BCALM.");
    parser.add("k", "K-mer length (must be <= " + std::to_string(constants::max_k) + ").");
    parser.add("m", "Minimizer length (must be < k and <= 32).");

    /* optional arguments */
    parser.add("seed", "Seed for minimizer computation (default is 42).", "-s", false);
    parser.add("threads", "Number of threads for pthash (default is 1).", "-t", false);
    parser.add("output_filename",
               "Output file name where the data structure will be serialized (no files generated "
               "by default).",
               "-o", false);
    parser.add("tmp_dirname",
               "Temporary directory used for construction in external memory (default is current "
               "directory '" +
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
               true);
    parser.add("verbose", "Verbose output during construction (disabled by default).", "--verbose",
               true);

    if (!parser.parse()) throw ParseError();
    return parser;
}

cmd_line_parser::parser get_query_parser(int argc, char* argv[]) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("mphf", "lphash minimal perfect hash function saved on disk\n");
    parser.add("input_filename",
               "Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not:\n");
    if (!parser.parse()) throw ParseError();
    return parser;
}

}  // namespace lphash
