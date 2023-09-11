#include <iostream>

#include "../include/constants.hpp"
#include "ptbb.hpp"

#include <argparse/argparse.hpp>

using namespace lphash;

int main(int argc, char* argv[]) {
    // cmd_line_parser::parser parser(argc, argv);
    argparse::ArgumentParser parser(argv[0]);
    parser.add_argument("-i", "--input-filename")
        .help("Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not:\n"
               "\t- without duplicate nor invalid kmers\n"
               "\t- one DNA sequence per line.\n"
               "\tFor example, it could be the de Bruijn graph topology output by BCALM.")
        .required();
    parser.add_argument("-k")
        .help("K-mer length (must be <= " + std::to_string(constants::max_k) + ").")
        .scan<'u', uint64_t>()
        .required();
    parser.add_argument("-p", "--pthash-filename")
        .help("Output file name where the pthash mphf will be serialized.");
    parser.add_argument("-b", "--bbhash-filename")
        .help("Output file name where the BBHash mphf will be serialized.");
    parser.add_argument("-a", "--alpha")
        .help("The table load factor. It must be a quantity > 0 and <= 1. (default is " + std::to_string(0.94) + ").")
        .scan<'g', double>()
        .default_value(double(0.94));
    parser.add_argument("-c")
        .help( "A (floating point) constant that trades construction speed for space effectiveness "
               "of minimal perfect hashing. "
               "A reasonable value lies between 3.0 and 10.0 (default is " +
                   std::to_string(constants::c) + ").")
        .scan<'g', double>()
        .default_value(double(constants::c));
    parser.add_argument("-g", "--gamma")
        .help("Load factor for BBHash (default is 1)")
        .scan<'g', double>()
        .default_value(double(1));
    parser.add_argument("-d", "--tmp-dirname")
        .help("Temporary directory used for construction in external memory. Default is directory '" + constants::default_tmp_dirname + "'.")
        .default_value(constants::default_tmp_dirname);
    parser.add_argument("-t", "--threads")
        .help("Number of threads (default is 1).")
        .scan<'u', std::size_t>()
        .default_value(std::size_t(1));
    parser.add_argument("-v", "--verbose")
        .help("Verbose output during construction.")
        .implicit_value(true)
        .default_value(false);
    parser.add_argument("-c", "--check")
        .help("Check output")
        .implicit_value(true)
        .default_value(false);
    try {
        parser.parse_args(argc, argv);
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    auto input_filename = parser.get<std::string>("-i");
    auto k = parser.get<uint64_t>("-k");

    gzFile fp = NULL;
    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
        std::cerr << "Unable to open the input file " << input_filename << "\n";
        return 2;
    }
    uint64_t total_kmers = 0;
    kseq_t* seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        std::string contig = std::string(seq->seq.s);
        uint64_t nbases_since_last_break = 0;
        for (uint64_t i = 0; i < contig.size(); ++i) {
            uint8_t c = constants::seq_nt4_table[static_cast<uint8_t>(contig[i])];
            if (c < 4) {
                ++nbases_since_last_break;
                if (nbases_since_last_break >= k) ++total_kmers;
            } else {
                nbases_since_last_break = 0;
            }
        }
    }
    if (seq) kseq_destroy(seq);
    gzclose(fp);
    std::cout << input_filename << "," << k << "," << total_kmers;

    ptbb::ptbb_file_itr kmer_itr(input_filename, k);
    ptbb::ptbb_file_itr kmer_end;
    uint64_t check_total_kmers = 0;
    for (; kmer_itr != kmer_end; ++kmer_itr) ++check_total_kmers;
    assert(total_kmers == check_total_kmers);

    auto num_threads = parser.get<std::size_t>("--threads");

    pthash::build_configuration pt_config;
    pt_config.num_threads = num_threads;

    bool check = parser.get<bool>("--check");
    pt_config.verbose_output = parser.get<bool>("--verbose");
    
    if (parser.is_used("-p")) {
        std::string pthash_filename = parser.get<std::string>("--pthash-filename");
        pt_config.minimal_output = true;
        pt_config.seed = constants::default_pthash_seed;
        pt_config.c = parser.get<double>("-c");
        pt_config.alpha = parser.get<double>("--alpha");
        pt_config.ram = 8 * essentials::GB;
        pt_config.tmp_dir = parser.get<std::string>("--tmp-dirname");
        if (pt_config.tmp_dir != constants::default_tmp_dirname) essentials::create_directory(pt_config.tmp_dir);
        ptbb::pthash_mphf_t pthash_mphf;
        {
            ptbb::ptbb_file_itr kmer_itr(input_filename, k);
            pthash_mphf.build_in_external_memory(kmer_itr, total_kmers, pt_config);
        }
        essentials::save(pthash_mphf, pthash_filename.c_str());
        assert(total_kmers == pthash_mphf.num_keys());
        std::cout << "," << pthash_mphf.num_bits() << ","
                  << static_cast<double>(pthash_mphf.num_bits()) / pthash_mphf.num_keys();

        if (check) {
            std::cerr << "Checking PTHash...";
            pthash::bit_vector_builder population(total_kmers);
            uint64_t check_total_kmers = 0;
            {
                ptbb::ptbb_file_itr kmer_itr(input_filename, k);
                ptbb::ptbb_file_itr kmer_end;
                for (; kmer_itr != kmer_end; ++kmer_itr) {
                    auto idx = pthash_mphf(*kmer_itr);
                    if (idx >= total_kmers) {
                        std::cerr << "[Error] out of bounds" << std::endl;
                        return 2;
                    } else if (population.get(idx)) {
                        std::cerr << "[Error] collision" << std::endl;
                        return 2;
                    } else
                        population.set(idx);
                    ++check_total_kmers;
                }
            }
            assert(total_kmers == check_total_kmers);
            for (uint64_t i = 0; i < total_kmers; ++i) {
                if (!population.get(i)) {
                    std::cerr << "[Error] hash is not perfect" << std::endl;
                    return 2;
                }
            }
            std::cerr << "EVERYTHING OK\n";
        }
    } else {
        std::cout << ",,";
    }

    if (parser.is_used("--bbhash-filename")) {
        std::string bbhash_filename = parser.get<std::string>("--bbhash-filename");
        double gammaFactor = parser.get<double>("--gamma");
        if (gammaFactor < 1.0) throw std::runtime_error("BBHash gamma factor < 1");

        std::vector<kmer_t> keys;
        ptbb::ptbb_file_itr boo_itr_begin(input_filename, k);
        ptbb::ptbb_file_itr boo_itr_end;
        for (; boo_itr_begin != boo_itr_end; ++boo_itr_begin) keys.push_back(*boo_itr_begin);
        auto data_iterator = boomphf::range(keys.begin(), keys.end());
        ptbb::bbhash_mphf_t bbhash_mphf(total_kmers, data_iterator, num_threads, gammaFactor, true,
                                        pt_config.verbose_output, 0);
        assert(keys.size() == total_kmers);
        keys.reserve(0);
        std::cout << "," << bbhash_mphf.totalBitSize() << ","
                  << static_cast<double>(bbhash_mphf.totalBitSize()) / total_kmers;

        if (check) {
            std::cerr << "Checking BBHash...";
            pthash::bit_vector_builder population(total_kmers);
            uint64_t check_total_kmers = 0;
            {
                ptbb::ptbb_file_itr kmer_itr(input_filename, k);
                ptbb::ptbb_file_itr kmer_end;
                for (; kmer_itr != kmer_end; ++kmer_itr) {
                    auto idx = bbhash_mphf.lookup(*kmer_itr);
                    if (idx >= total_kmers) {
                        std::cerr << "[Error] out of bounds: " << idx << std::endl;
                        return 2;
                    } else if (population.get(idx)) {
                        std::cerr << "[Error] collision" << std::endl;
                        return 2;
                    } else
                        population.set(idx);
                    ++check_total_kmers;
                }
            }
            assert(total_kmers == check_total_kmers);
            for (uint64_t i = 0; i < total_kmers; ++i) {
                if (!population.get(i)) {
                    std::cerr << "[Error] hash is not perfect" << std::endl;
                    return 2;
                }
            }
            std::cerr << "EVERYTHING OK\n";
        }

        {
            std::ofstream bbh_strm(bbhash_filename, std::ios::binary);
            bbhash_mphf.save(bbh_strm);
        }
    } else {
        std::cout << ",,";
    }

    std::cout << std::endl;
}