extern "C" {
#include "../include/kseq.h"
}
#include <zlib.h>
#include <string>

#include "../include/constants.hpp"
#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "../include/ptbb_file_itr.hpp"
#include "../include/BooPHF.hpp"

KSEQ_INIT(gzFile, gzread)

using namespace lphash;

int main(int argc, char* argv[]) {
    gzFile fp;
    kseq_t* seq;
    int c;
    std::size_t total_kmers;
    uint64_t k;
    std::string input_filename;
    bool check;
    
    cmd_line_parser::parser parser(argc, argv);
    parser.add("input_filename",
               "Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not:\n"
               "\t- without duplicate nor invalid kmers\n"
               "\t- one DNA sequence per line.\n"
               "\tFor example, it could be the de Bruijn graph topology output by BCALM.");
    parser.add("k", "K-mer length (must be <= " + std::to_string(constants::max_k) + ").");
    parser.add("pthash_filename", 
               "Output file name where the pthash mphf will be serialized.",
               "-p",
               false);
    parser.add("bbhash_filename", 
               "Output file name where the BBHash mphf will be serialized.",
               "-b",
               false);
    parser.add("alpha",
               "The table load factor. It must be a quantity > 0 and <= 1. (default is " + std::to_string(0.94) + ").",
               "-a",
               false);
    parser.add("gamma",
               "Load factor for BBHash (default is 1)",
               "-g",
               false);
    parser.add("tmp_dirname",
               "Temporary directory used for construction in external memory. Default is directory '" + constants::default_tmp_dirname + "'.",
               "-d", 
               false);
    parser.add("threads", 
               "Number of threads for pthash (default is 1).", 
               "-t", 
               false);
    parser.add("verbose", "Verbose output during construction.", "--verbose", true);
    parser.add("check", "Check output", "--check", true);
    if (!parser.parse()) return 1;

    input_filename = parser.get<std::string>("input_filename");
    k = parser.get<uint64_t>("k");

    fp = NULL;
    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
        std::cerr << "Unable to open the input file " << input_filename << "\n";
        return 2;
    }
    total_kmers = 0;
    seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        std::string contig = std::string(seq->seq.s);  // we lose a little bit of efficiency here
        std::size_t nbases_since_last_break = 0;
        for (std::size_t i = 0; i < contig.size(); ++i) {
            c = constants::seq_nt4_table[static_cast<uint8_t>(contig[i])];
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

    other::ptbb_file_itr kmer_itr(input_filename, k);
    // kseq_t* itr_guts = reinterpret_cast<kseq_t*>(kmer_itr.memory_management());
    other::ptbb_file_itr kmer_end;
    std::size_t check_total_kmers = 0;
    for (; kmer_itr != kmer_end; ++kmer_itr) {
        ++check_total_kmers;
    }
    // kseq_destroy(itr_guts);
    std::cerr << "\n" << total_kmers << " == " << check_total_kmers << std::endl;
    assert(total_kmers == check_total_kmers);

    pthash::build_configuration pt_config;
    if (parser.parsed("threads")) pt_config.num_threads = parser.get<uint32_t>("threads");
    else pt_config.num_threads = 1;
    if (parser.parsed("check")) check = true;
    else check = false;
    if (parser.parsed("pthash_filename")) {
        std::string pthash_filename = parser.get<std::string>("pthash_filename");
        pt_config.minimal_output = true;
        pt_config.seed = constants::seed;
        pt_config.c = constants::c;
        pt_config.alpha = 0.94;
        if (parser.parsed("verbose")) pt_config.verbose_output = parser.get<bool>("verbose");
        else pt_config.verbose_output = false;
        if (parser.parsed("tmp_dirname")) {
            pt_config.tmp_dir = parser.get<std::string>("tmp_dirname");
            essentials::create_directory(pt_config.tmp_dir);
        }
        pthash_mphf_t kmer_order;
        {
            other::ptbb_file_itr kmer_itr(input_filename, k);
            // kseq_t* itr_guts = reinterpret_cast<kseq_t*>(kmer_itr.memory_management());
            kmer_order.build_in_external_memory(kmer_itr, total_kmers, pt_config);
        }
        // kseq_destroy(itr_guts);
        essentials::save(kmer_order, pthash_filename.c_str());
        
        assert(total_kmers == kmer_order.num_keys());
        std::cout << "," << kmer_order.num_bits() << "," << static_cast<double>(kmer_order.num_bits()) / kmer_order.num_keys();

        if (check) {
            std::cerr << "Checking PTHash" << std::endl;
            pthash::bit_vector_builder population(total_kmers);
            std::size_t check_total_kmers = 0;
            {
                other::ptbb_file_itr kmer_itr(input_filename, k);
                // kseq_t* itr_guts = reinterpret_cast<kseq_t*>(kmer_itr.memory_management());
                other::ptbb_file_itr kmer_end;
                for (; kmer_itr != kmer_end; ++kmer_itr) {
                    auto idx = kmer_order(*kmer_itr);
                    if (idx >= total_kmers) {
                        std::cerr << "[Error] out of bounds" << std::endl;
                        return 2;
                    } else if (population.get(idx)) {
                        std::cerr << "[Error] collision" << std::endl;
                        return 2;
                    } else population.set(idx);
                    ++check_total_kmers;
                }
            } // kseq_destroy(itr_guts);
            assert(total_kmers == check_total_kmers);
            for (std::size_t i = 0; i < total_kmers; ++i) {
                if (!population.get(i)) {
                    std::cerr << "[Error] hash is not perfect" << std::endl;
                    return 2;
                }
            }
        }
    }
    if (parser.parsed("bbhash_filename")) {
        std::string bbhash_filename = parser.get<std::string>("bbhash_filename");
        double gammaFactor; // lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query
        if (parser.parsed("gamma")) gammaFactor = parser.get<double>("gamma");
        else gammaFactor = 1.0;
        if (gammaFactor < 1.0) throw std::runtime_error("BBHash gamma factor < 1");
        other::ptbb_file_itr boo_itr_begin(input_filename, k);
        other::ptbb_file_itr boo_itr_end;
        auto data_iterator = boomphf::range(boo_itr_begin, boo_itr_end);
        boomphf::mphf<kmer_t, other::BBHasher<kmer_t>> bphf(total_kmers, data_iterator, pt_config.num_threads, gammaFactor);
        
        std::cout << "," << bphf.totalBitSize() << "," << static_cast<double>(bphf.totalBitSize()) / total_kmers;

        if (check) {
            std::cerr << "Checking BBHash" << std::endl;
            
            pthash::bit_vector_builder population(total_kmers);
            std::size_t check_total_kmers = 0;
            // other::BBHasher<kmer_t> bbhasher;
            {
                other::ptbb_file_itr kmer_itr(input_filename, k);
                // kseq_t* itr_guts = reinterpret_cast<kseq_t*>(kmer_itr.memory_management());
                other::ptbb_file_itr kmer_end;
                for (; kmer_itr != kmer_end; ++kmer_itr) {
                    // std::cerr << "Hash = " << bbhasher(*kmer_itr) << std::endl;
                    auto idx = bphf.lookup(*kmer_itr);
                    if (idx >= total_kmers) {
                        std::cerr << "[Error] out of bounds: " << idx << std::endl;
                        return 2;
                    } else if (population.get(idx)) {
                        std::cerr << "[Error] collision" << std::endl;
                        return 2;
                    } else population.set(idx);
                    ++check_total_kmers;
                }
            }
            // kseq_destroy(itr_guts);
            assert(total_kmers == check_total_kmers);
            for (std::size_t i = 0; i < total_kmers; ++i) {
                if (!population.get(i)) {
                    std::cerr << "[Error] hash is not perfect" << std::endl;
                    return 2;
                }
            }
        }

        {//
            std::ofstream bbh_strm(bbhash_filename, std::ios::binary);
            bphf.save(bbh_strm);
        }//
    }
   std::cout << "\n";
}