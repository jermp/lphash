#pragma once

#include <string>

namespace lphash {

#include "../compile_constants.tpd"

namespace constants {

/* max *odd* size that can be packed into the given k-mer type */
constexpr uint64_t max_k = (sizeof(kmer_t) * 8) / 2 - 1;
constexpr uint64_t default_pthash_seed = 1;
constexpr uint64_t default_seed = 42;
constexpr uint64_t default_num_threads = 1;
constexpr double c = 3.0;  // for PTHash
static const std::string default_tmp_dirname(".");

}

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

}