#pragma once

#include <exception>
#include <string>

namespace lphash {

class ParseError : public std::exception {
public:
    const char* what() const noexcept { return "Unable to parse the arguments\n"; }
};

class OptionError : public std::exception {
public:
    OptionError(const char* message) : msg(message) {}
    OptionError(std::string const& message) : msg(message) {}
    const char* what() const noexcept { return msg.c_str(); }
private:
    std::string msg;
};

struct configuration {
    configuration();
    std::string input_filename;
    std::string output_filename;
    uint64_t k;
    uint64_t m;
    uint64_t mm_seed;
    double c;
    uint64_t num_threads;
    uint64_t max_memory;
    std::string tmp_dirname;
    bool check, verbose;
};

}