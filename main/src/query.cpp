#include <argparse/argparse.hpp>

namespace lphash {

argparse::ArgumentParser get_parser_query() {
    argparse::ArgumentParser parser("query");
    parser.add_description("query a LP-MPHF");
    parser.add_argument("-i", "--input-mphf")
        .help("LP-Hash MPHF saved on disk")
        .required();
    parser.add_argument("-q", "--query-filename")
        .help("Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not")
        .required();
    return parser;
}

}  // namespace lphash