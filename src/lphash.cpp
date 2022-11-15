#include <iostream>
#include "../include/build.hpp"
#include "../include/query.hpp"

using namespace lphash;

int help(char* arg0) {
    std::cerr << "LP-Hash: Locality Preserving Minimal Perfect Hahsing of k-mers\n\n";
    std::cerr << "Usage: " << arg0 << " <tool> ...\n\n";
    std::cerr << "Available tools:\n";
    std::cerr << "  build-p      build a partitioned LP-MPHF\n";
    std::cerr << "  build-u      build an unpartitioned LP-MPHF\n";
    std::cerr << "  query-p      query a partitioned LP-MPHF\n";
    std::cerr << "  query-u      query an unpartitioned LP-MPHF\n";
    return 1;
}

int main(int argc, char* argv[]) {
    if (argc < 2) return help(argv[0]);
    auto tool = std::string(argv[1]);
    if (tool == "build-p") {
        return build_partitioned_main(argc - 1, argv + 1);
    } else if (tool == "build-u") {
        return build_unpartitioned_main(argc - 1, argv + 1);
    } else if (tool == "query-p") {
        return query_partitioned_main(argc - 1, argv + 1);
    } else if (tool == "query-u") {
        return query_unpartitioned_main(argc - 1, argv + 1);
    }
    std::cerr << "Unsupported tool '" << tool << "'." << std::endl;
    return help(argv[0]);
}