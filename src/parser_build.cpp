#include "../include/constants.hpp"
#include "../include/parser_build.hpp"

namespace lphash {

cmd_line_parser::parser get_build_parser(int argc, char* argv[]) {
    cmd_line_parser::parser parser(argc, argv);

    /* mandatory arguments */
    parser.add("input_filename",
               "Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not:\n"
               "\t- without duplicate nor invalid kmers\n"
               "\t- one DNA sequence per line.\n"
               "\tFor example, it could be the de Bruijn graph topology output by BCALM.");
    parser.add("k", "K-mer length (must be <= " + std::to_string(constants::max_k) + ").");
    parser.add("m", "Minimizer length (must be < k).");

    /* optional arguments */
    parser.add("seed",
               "Seed for minimizer computation (default is 42).", 
               "-s",
               false);
    // parser.add("c",
    //            "A (floating point) constant that trades construction speed for space effectiveness of minimal perfect hashing. "
    //            "A reasonable value lies between 3.0 and 10.0 (default is " + std::to_string(constants::c) + ").",
    //            "-c", 
    //            false);
    // parser.add("a", 
    //            "(default is 0.94 ).", 
    //            "-a", 
    //            false);
    parser.add("output_filename", 
               "Output file name where the data structure will be serialized.",
               "-o", 
               false);
    parser.add("results_filename",
               "CSV file to append to",
               "-r",
               false);
    parser.add("tmp_dirname",
               "Temporary directory used for construction in external memory. Default is directory '" + constants::default_tmp_dirname + "'.",
               "-d", 
               false);
    parser.add("canonical_parsing",
               "Canonical parsing of k-mers. "
               "This option changes the parsing and results in a trade-off between index space and lookup time.",
               "--canonical-parsing", 
               true);
    parser.add("check", "Check correctness after construction.", "--check", true);
    // parser.add("bench", "Run benchmark after construction.", "--bench", true);
    parser.add("verbose", "Verbose output during construction.", "--verbose", true);

    if (!parser.parse()) throw ParseError();
    return parser;
}

}
