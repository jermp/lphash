#ifndef PARSER_BUILD_HPP
#define PARSER_BUILD_HPP

#include <exception>
#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"

namespace lphash {

class ParseError : public std::exception {
    public:
        const char * what() const noexcept { return "Unable to parse the arguments\n"; }
};

cmd_line_parser::parser get_build_parser(int argc, char* argv[]);

}

#endif // PARSER_BUILD_HPP