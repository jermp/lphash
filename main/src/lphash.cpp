#include <zlib.h>
extern "C" {
#include "../../external/kseq.h"
}
#include <iostream>
#include "../../lib/include/partitioned_mphf.hpp"
#include "../../lib/include/unpartitioned_mphf.hpp"

KSEQ_INIT(gzFile, gzread)
#include "../include/build.hpp"
#include "../include/query.hpp"

using namespace lphash;

argparse::ArgumentParser get_parser_partitioned();
argparse::ArgumentParser get_parser_unpartitioned();

template<class MPHF>
int dispatch(const argparse::ArgumentParser& parser) {
    if (parser.is_subcommand_used("build")) return build_main<MPHF>(parser);
    else if (parser.is_subcommand_used("query")) return query_main<MPHF>(parser);
    else throw std::runtime_error("This should never happen");
}

int main(int argc, char* argv[])
{
    auto build_parser = get_parser_build();
    auto query_parser = get_parser_query();
    auto pparser = argparse::ArgumentParser("partitioned");
    pparser.add_description("Build a partitioned LP-Hash");
    auto uparser = argparse::ArgumentParser("unpartitioned");
    uparser.add_description("Build an unpartitioned LP-Hash");

    pparser.add_subparser(build_parser);
    pparser.add_subparser(query_parser);
    uparser.add_subparser(build_parser);
    uparser.add_subparser(query_parser);

    argparse::ArgumentParser program(argv[0]);
    program.add_description("LP-Hash: (L)ocality (P)reserving Minimal Perfect (Hash)ing of k-mers");
    
    program.add_subparser(pparser);
    program.add_subparser(uparser);
    try {
        program.parse_args(argc, argv);
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << program;
        return 1;
    }

    if (program.is_subcommand_used(pparser)) return dispatch<mphf::partitioned>(pparser);
    else if (program.is_subcommand_used(uparser)) return dispatch<mphf::unpartitioned>(uparser);
    else std::cerr << program << std::endl;
    return 0;
}