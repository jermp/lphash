#include <iostream>
#include "../include/build.hpp"
#include "../include/query.hpp"

using namespace lphash;

void print_commands() {
    std::cerr << "Available commands:\n";
    std::cerr << "\tpartitioned\tbuild a partitioned lp-mphf\n";
    std::cerr << "\tunpartitioned\tbuild an unpartitioned lp-mphf\n";
    std::cerr << "\tqp\tquery a partitioned lp-mphf\n";
    std::cerr << "\tqu\tquery an unpartitioned lp-mphf\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_commands();
        return 1;
    }
    auto cmd = std::string(argv[1]);
    if (cmd == "partitioned") {
        return build_partitioned_main(argc - 1, argv + 1);
    } else if (cmd == "unpartitioned") {
        return build_unpartitioned_main(argc - 1, argv + 1);
    } else if (cmd == "qp") {
        return query_partitioned_main(argc - 1, argv + 1);
    } else if (cmd == "qu") {
        return query_unpartitioned_main(argc - 1, argv + 1);
    } else {
        std::cerr << "Unsupported command '" << cmd << "'." << std::endl;
        print_commands();
        return 1;
    }
}