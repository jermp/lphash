extern "C" {
    #include "../include/kseq.h"
}

#include <zlib.h>
#include "../include/constants.hpp"
#include "../include/quartet_wtree.hpp"
#include "../external/pthash/external/cmd_line_parser/include/parser.hpp"
#include "minimizer.hpp"

#include "../include/prettyprint.hpp"

using namespace lphash;

typedef pthash::murmurhash2_64 base_hasher_type;
typedef pthash::single_phf<base_hasher_type,               // base hasher
                           pthash::dictionary_dictionary,  // encoder type
                           true                            // minimal output
                           >
    pthash_mphf_type;

KSEQ_INIT(gzFile, gzread)

class vector_mmp_to_pthash_itr_adapter : std::forward_iterator_tag {
    public:
        typedef mmp_t value_type;
        vector_mmp_to_pthash_itr_adapter(std::vector<mmp_t>::iterator begin, std::vector<mmp_t>::iterator end) : begin(begin), end(end), current(begin) {};
        // inline uint64_t minimizer() const {return current->itself;};
        inline uint64_t operator*() const {return current->itself;};
        inline void operator++() {
            uint64_t prev_mm = current->itself;
            while(current != end && current->itself == prev_mm) {++current;}
        };
    private:
        std::vector<mmp_t>::iterator begin, end, current;
};

int main (int argc, char* argv[]) 
{
    gzFile fp;
    kseq_t *seq;
    cmd_line_parser::parser parser(argc, argv);

    /* mandatory arguments */
    parser.add("input_filename",
               "Must be a FASTA file (.fa/fasta extension) compressed with gzip (.gz) or not:\n"
               "\t- without duplicate nor invalid kmers\n"
               "\t- one DNA sequence per line.\n"
               "\tFor example, it could be the de Bruijn graph topology output by BCALM.");
    parser.add("k", "K-mer length (must be <= " + std::to_string(constants::max_k) + ").");
    parser.add("m", "Minimizer length (must be < k).");

    /* optional arguments */
    parser.add("seed",
               "Seed for construction (default is " + std::to_string(constants::seed) + ").", "-s",
               false);
    parser.add("l",
               "A (integer) constant that controls the space/time trade-off of the dictionary. "
               "A reasonable values lies between 2 and 12 (default is " +
                   std::to_string(constants::min_l) + ").",
               "-l", false);
    parser.add("c",
               "A (floating point) constant that trades construction speed for space effectiveness "
               "of minimal perfect hashing. "
               "A reasonable value lies between 3.0 and 10.0 (default is " +
                   std::to_string(constants::c) + ").",
               "-c", false);
    parser.add("output_filename", "Output file name where the data structure will be serialized.",
               "-o", false);
    parser.add(
        "tmp_dirname",
        "Temporary directory used for construction in external memory. Default is directory '" +
            constants::default_tmp_dirname + "'.",
        "-d", false);
    parser.add("canonical_parsing",
               "Canonical parsing of k-mers. This option changes the parsing and results in a "
               "trade-off between index space and lookup time.",
               "--canonical-parsing", true);
    parser.add("check", "Check correctness after construction.", "--check", true);
    parser.add("bench", "Run benchmark after construction.", "--bench", true);
    parser.add("verbose", "Verbose output during construction.", "--verbose", true);

    if (!parser.parse()) return 1;

    auto input_filename = parser.get<std::string>("input_filename");
    fp = NULL;
    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
        std::cerr << "Unable to open the input file " << input_filename << "\n";
        return 2;
    }

    if(fp == NULL) {// read file from stdin (not used as of now)
        if ((fp = gzdopen(fileno(stdin), "r")) == NULL) {
            fprintf(stderr, "Unable to use stdin as input\n");
            return 2;
        }
    }

    auto k = parser.get<uint32_t>("k");
    auto m = parser.get<uint32_t>("m");
    auto seed = parser.get<uint64_t>("seed");
    std::size_t total_kmers = 0;
    std::vector<mmp_t> minimizers;
    seq = kseq_init(fp);
    
    /*part 1: read sequences, add to each minimizer its position in the first k-mer of the super-k-mer and the length of the super-k-mer*/
    std::cerr << "Part 1: file reading and info gathering" << std::endl;
    while(kseq_read(seq) >= 0) {
        std::string contig = std::string(seq->seq.s); //we loose a little bit of efficiency here
        // debug::print_hashes(contig, m, seed);
        auto n = minimizer::from_string<murmurhash2_64>(contig, k, m, seed, false, minimizers);//not canonical minimizers for now
        total_kmers += n;
        std::cout << "read " << n << " k-mers (contig length = " << contig.length() << ")\n";
        // debug::compute_minimizers_naive<murmurhash2_64>(contig, k, m, seed);
    }
    if (seq) kseq_destroy(seq);
    // std::cerr << std::endl;
    // for (auto& mm : minimizers) std::cerr << mm << std::endl;
    
    /*part 2: build MPHF from minimizers*/
    std::cerr << "Part 2: build minimizer MPHF" << std::endl;
    auto mm_compare = [](mmp_t const& a, mmp_t const& b) {return a.itself < b.itself;};
    std::sort(minimizers.begin(), minimizers.end(), mm_compare);
    std::cerr << "--- : minimizers are now sorted by their value" << std::endl;
    std::size_t n_distinct_minimizers = 0;
    for(auto it = minimizers.begin(), prev = minimizers.end(); it != minimizers.end(); ++it) {
        if (prev->itself != it->itself) {
            ++n_distinct_minimizers;
            prev = it;
        }
    }
    std::cerr << "--- : [Warning] get number of distinct minimizer -> rework pthash interface for this" << std::endl;
    pthash_mphf_type mm_mphf;
    pthash::build_configuration mphf_config;
    mphf_config.c = 6.0;
    mphf_config.alpha = 0.94;
    mphf_config.seed = 42;  // my favourite seed
    mphf_config.minimal_output = true;
    mphf_config.verbose_output = true;
    mphf_config.num_threads = 1;
    uint64_t num_threads = std::thread::hardware_concurrency() >= 8 ? 8 : 1;
    // if (minimizers.size() >= num_threads) mphf_config.num_threads = num_threads;
    mphf_config.ram = 2 * essentials::GB;
    if (parser.parsed("tmp_dirname")) {
        mphf_config.tmp_dir = parser.get<std::string>("tmp_dirname");
        essentials::create_directory(mphf_config.tmp_dir);
    }
    std::cerr << "--- : control 1" << std::endl;
    auto begin = vector_mmp_to_pthash_itr_adapter(minimizers.begin(), minimizers.end());
    std::cerr << "--- : number of minimizers: " << minimizers.size() << ", of which distinct: " << n_distinct_minimizers << std::endl;
    mm_mphf.build_in_external_memory(begin, n_distinct_minimizers, mphf_config);
    std::cerr << "--- : minimizers are now sorted by MPFH" << std::endl;

    /*part 3: sort vector based on MPHF order*/
    std::cerr << "Part 3: sort minimizers by MPHF" << std::endl;
    auto mphf_compare = [&mm_mphf](mmp_t const& a, mmp_t const& b) {return mm_mphf(a.itself) < mm_mphf(b.itself);};
    std::sort(minimizers.begin(), minimizers.end(), mphf_compare);

    /*part 4: build wavelet tree and <positions/size> vectors*/
    std::cerr << "Part 4: build wavelet tree and offsets" << std::endl;
    quartet_wtree_builder wtb(minimizers.size());
    std::vector<uint64_t> colliding;
    uint64_t n_maximal = 0;
    std::vector<uint64_t> right_or_collision_sizes;
    std::vector<uint64_t> left_positions;
    std::vector<uint64_t> none_positions, none_sizes;
    for (std::size_t i = 0; i < minimizers.size(); ++i) {
        // check if minimizer there are multiple minimizers
        std::size_t j;
        for (j = i; j < minimizers.size() && minimizers[j].itself == minimizers[i].itself; ++j) {}
        if ((j-i) == 1) {
            if (minimizers[i].p1 == k-m+1) {
                if (minimizers[i].size == k-m+1) {
                    wtb.push_back(MAXIMAL);
                    ++n_maximal;
                } else {
                    wtb.push_back(RIGHT_OR_COLLISION);
                    right_or_collision_sizes.push_back(minimizers[i].size);
                }
            } else {
                if (minimizers[i].p1 == minimizers[i].size - 1) {
                    wtb.push_back(LEFT);
                    left_positions.push_back(minimizers[i].p1);
                } else {
                    wtb.push_back(NONE);
                    none_positions.push_back(minimizers[i].p1);
                    none_sizes.push_back(minimizers[i].size);
                }
            }
        } else { // collision
            wtb.push_back(RIGHT_OR_COLLISION);
            colliding.push_back(minimizers[i].itself);
            right_or_collision_sizes.push_back(0);
        }
    }
    assert(none_positions.size() == none_sizes.size());
    std::cerr << "Number of Maximal minimizers: " << n_maximal << "\n";
    std::cerr << "Number of Leftmax minimizers: " << left_positions.size() << "\n";
    std::cerr << "Number of Rightmax minimizers: " << right_or_collision_sizes.size() << "\n";
    std::cerr << "Number of Uncategorized minimizers: " << none_positions.size() << "\n";
    std::cerr << "Number of Ambiguous minimizers: " << colliding.size() << "\n";

    /*part 5: Elias-Fano*/

    /*part 6: build fallback mphf*/

    /*part 7: save everything*/

    /*part 8: check for correctness (query + check minimality)*/
}