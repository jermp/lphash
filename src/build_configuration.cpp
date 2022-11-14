#include "../include/build_configuration.hpp"
#include "../include/constants.hpp"

namespace lphash {

configuration::configuration() 
    : input_filename(""), 
      output_filename(""),
      k(0), 
      m(0), 
      mm_seed(constants::default_seed), 
      c(constants::c), 
      num_threads(constants::default_num_threads),
      max_memory(8),
      tmp_dirname(constants::default_tmp_dirname),
      check(false), verbose(false) {}

} // namespace lphash