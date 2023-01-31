#include <iostream>

#include "../../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "ptbb.hpp"

using namespace lphash;

int main(int argc, char* argv[]) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("input_filename",
               "Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not:\n"
               "\t- without duplicate nor invalid kmers\n"
               "\t- one DNA sequence per line.\n"
               "\tFor example, it could be the de Bruijn graph topology output by BCALM.",
               "-i", true);
    parser.add("k", "K-mer length (must be <= " + std::to_string(constants::max_k) + ").", "-k",
               true);
    parser.add("pthash_filename", "Output file name where the pthash mphf will be serialized.",
               "-p", false);
    parser.add("bbhash_filename", "Output file name where the BBHash mphf will be serialized.",
               "-b", false);
    parser.add("alpha",
               "The table load factor. It must be a quantity > 0 and <= 1. (default is " +
                   std::to_string(0.94) + ").",
               "-a", false);
    parser.add("c",
               "A (floating point) constant that trades construction speed for space effectiveness "
               "of minimal perfect hashing. "
               "A reasonable value lies between 3.0 and 10.0 (default is " +
                   std::to_string(constants::c) + ").",
               "-c", false);
    parser.add("gamma", "Load factor for BBHash (default is 1)", "-g", false);
    parser.add(
        "tmp_dirname",
        "Temporary directory used for construction in external memory. Default is directory '" +
            constants::default_tmp_dirname + "'.",
        "-d", false);
    parser.add("threads", "Number of threads (default is 1).", "-t", false);
    parser.add("verbose", "Verbose output during construction.", "--verbose", false, true);
    parser.add("check", "Check output", "--check", false, true);
    if (!parser.parse()) return 1;

    auto input_filename = parser.get<std::string>("input_filename");
    auto k = parser.get<uint64_t>("k");

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
    for (; kmer_itr != kmer_end; ++kmer_itr) { ++check_total_kmers; }
    assert(total_kmers == check_total_kmers);

    uint32_t num_threads = 1;
    if (parser.parsed("threads")) num_threads = parser.get<uint32_t>("threads");

    pthash::build_configuration pt_config;
    pt_config.num_threads = num_threads;

    bool check = parser.get<bool>("check");
    if (parser.parsed("verbose")) {
        pt_config.verbose_output = parser.get<bool>("verbose");
    } else {
        pt_config.verbose_output = false;
    }

    if (parser.parsed("pthash_filename")) {
        std::string pthash_filename = parser.get<std::string>("pthash_filename");
        pt_config.minimal_output = true;
        pt_config.seed = constants::default_pthash_seed;
        pt_config.c = (parser.parsed("c")) ? parser.get<double>("c") : constants::c;
        pt_config.alpha = 0.94;
        pt_config.ram = 8 * essentials::GB;
        if (parser.parsed("tmp_dirname")) {
            pt_config.tmp_dir = parser.get<std::string>("tmp_dirname");
            essentials::create_directory(pt_config.tmp_dir);
        }
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

    if (parser.parsed("bbhash_filename")) {
        std::string bbhash_filename = parser.get<std::string>("bbhash_filename");
        double gammaFactor;
        if (parser.parsed("gamma")) {
            gammaFactor = parser.get<double>("gamma");
        } else {
            gammaFactor = 1.0;
        }
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