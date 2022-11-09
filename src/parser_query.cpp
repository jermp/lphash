#include "../include/constants.hpp"
#include "../include/parser_build.hpp"

namespace lphash {

cmd_line_parser::parser get_query_parser(int argc, char* argv[]) {
    cmd_line_parser::parser parser(argc, argv);

    /* mandatory arguments */
    parser.add("mphf", "lphash minimal perfect hash function saved on disk\n");
    parser.add("input_filename",
               "Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not:\n");

    /* optional arguments */
    // parser.add("output_filename", 
    //            "Output file name with the result of the query (default is stdout)",
    //            "-o", 
    //            false);
    // parser.add("results_filename",
    //            "CSV file to append to",
    //            "-r",
    //            false);
    // parser.add("canonical_parsing",
    //            "Canonical parsing of k-mers. "
    //            "This option changes the parsing and results in a trade-off between index space and lookup time.",
    //            "--canonical-parsing", 
    //            true);
    if (!parser.parse()) throw ParseError();
    return parser;
}

}