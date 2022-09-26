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
    // parser.add("canonical_parsing",
    //            "Canonical parsing of k-mers. This option changes the parsing and results in a "
    //            "trade-off between index space and lookup time.",
    //            "--canonical-parsing", true);
    // parser.add("check", "Check correctness after construction.", "--check", true);
    // parser.add("bench", "Run benchmark after construction.", "--bench", true);
    // parser.add("verbose", "Verbose output during construction.", "--verbose", true);

    if (!parser.parse()) return 1;

    auto input_filename = parser.get<std::string>("input_filename");
    fp = NULL;
    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) {
        std::cerr << "Unable to open the input file " << input_filename << "\n";
        return 2;
    }

    auto k = parser.get<uint32_t>("k");
    auto m = parser.get<uint32_t>("m");
    auto seed = parser.get<uint64_t>("seed");
    std::size_t total_kmers = 0, total_kmers_check = 0;
    std::vector<mmp_t> minimizers;
    seq = kseq_init(fp);
    /*part 1: read sequences, add to each minimizer its position in the first k-mer of the super-k-mer and the length of the super-k-mer*/
    std::cout << "Part 1: file reading and info gathering\n";
    std::size_t old_size = 0;
    while(kseq_read(seq) >= 0) {
        std::string contig = std::string(seq->seq.s); //we lose a little bit of efficiency here
        auto n = minimizer::from_string<murmurhash2_64>(contig, k, m, seed, false, minimizers);//not canonical minimizers for now
        auto n_check = contig.length() - k + 1;
        total_kmers_check += n_check;
        // std::cerr << "read " << n << " k-mers (number of k-mers in contig = " << n_check << ")\n";
        // for (std::size_t i = old_size; i < minimizers.size(); ++i) std::cerr << minimizers[i].itself << "\n";
        // old_size = minimizers.size();
        // if (n != n_check) {
        //     std::cout << contig << "\n";
        // }
        // debug::compute_minimizers_naive<murmurhash2_64>(contig, k, m, seed);
        // debug::print_hashes(contig, m, seed);
        total_kmers += n;
    }
    if (seq) kseq_destroy(seq);
    assert(total_kmers == total_kmers_check); // valid iff contigs do not contain Ns
    // {
    //     std::ofstream mmfile("minimizers.txt");
    //     for (auto mm : minimizers) mmfile << mm << "\n";
    // }
    
    /*part 2: build MPHF from minimizers*/
    std::cout << "Part 2: build minimizer MPHF\n";
    auto mm_compare = [](mmp_t const& a, mmp_t const& b) {return a.itself < b.itself;};
    std::sort(minimizers.begin(), minimizers.end(), mm_compare);
    std::cout << "--- : minimizers are now sorted by their value\n";
    // std::set<uint64_t> mm_set;
    std::cout << "--- : [Warning] get number of distinct minimizer -> rework pthash interface for this\n";
    std::size_t n_distinct_minimizers = 0;
    for(auto it = minimizers.begin(), prev = minimizers.begin(); it != minimizers.end(); ++it) {
        // mm_set.insert(it->itself);
        if (prev->itself != it->itself) {
            ++n_distinct_minimizers;
            prev = it;
        }
    }
    if (minimizers.size()) ++n_distinct_minimizers;
    // std::cerr << n_distinct_minimizers << " vs " << mm_set.size() << "\n";
    // assert(n_distinct_minimizers == mm_set.size());
    pthash_mphf_type mm_mphf;
    pthash::build_configuration mphf_config;
    mphf_config.c = constants::c;
    mphf_config.seed = 42; // my favourite seed, different from the minimizer's seed.
    // mphf_config.seed = constants::seed;  
    if (parser.parsed("c")) mphf_config.c = parser.get<double>("c");
    mphf_config.alpha = 0.99;//0.94;
    mphf_config.minimal_output = true;
    mphf_config.verbose_output = false;
    mphf_config.num_threads = 1; // minimizers.size() >= std::thread::hardware_concurrency() ? std::thread::hardware_concurrency() : 1;
    mphf_config.ram = 2 * essentials::GB;
    // std::cerr << "n threads = " << mphf_config.num_threads << "\n";
    if (parser.parsed("tmp_dirname")) {
        mphf_config.tmp_dir = parser.get<std::string>("tmp_dirname");
        essentials::create_directory(mphf_config.tmp_dir);
    }
    auto begin = vector_mmp_to_pthash_itr_adapter(minimizers.begin(), minimizers.end());
    std::cout << "--- : number of minimizers: " << minimizers.size() << ", of which distinct: " << n_distinct_minimizers << "\n";
    mm_mphf.build_in_external_memory(begin, n_distinct_minimizers, mphf_config);

    /*part 3: sort vector based on MPHF order*/
    std::cout << "Part 3: sort minimizers by MPHF\n";
    auto mphf_compare = [&mm_mphf](mmp_t const& a, mmp_t const& b) {return mm_mphf(a.itself) < mm_mphf(b.itself);};
    std::sort(minimizers.begin(), minimizers.end(), mphf_compare);

    /*part 4: build wavelet tree and <positions/size> vectors*/
    std::cout << "Part 4: build wavelet tree and offsets\n";
    // pthash::bit_vector_builder test;
    // test.reserve(3);
    // test.push_back(true); test.push_back(false); test.push_back(true);
    // rs_bit_vector bv;
    // bv.build(&test, false);
    // std::cerr << bv[0] << bv[1] << bv[2] << std::endl;
    // std::cerr << bv.rank(0) << " " << bv.rank(1) << " " << bv.rank(2) << " " << bv.rank(3) << "\n";
    // std::cerr << bv.rank0(0) << " " << bv.rank0(1) << " " << bv.rank0(2) << " " << bv.rank0(3) << "\n";

    quartet_wtree_builder wtb(minimizers.size());
    std::vector<uint64_t> colliding_minimizers;
    uint64_t n_maximal = 0;
    std::vector<uint64_t> right_or_collision_sizes;
    std::vector<uint64_t> left_positions;
    std::vector<uint64_t> none_positions, none_sizes;
    for (std::size_t i = 0; i < minimizers.size(); ++i) {
        if (minimizers[i].itself != minimizers[i+1].itself) { // unique minimizer?
            if (minimizers[i].p1 == k-m) {
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
                    left_positions.push_back(minimizers[i].p1 + 1); // +1 because with p1 == 0 we have 1 k-mer in the prefix sum
                } else {
                    wtb.push_back(NONE);
                    none_positions.push_back(minimizers[i].p1);
                    none_sizes.push_back(minimizers[i].size);
                }
            }
        } else { // collision
            wtb.push_back(RIGHT_OR_COLLISION);
            colliding_minimizers.push_back(minimizers[i].itself);
            right_or_collision_sizes.push_back(0);
            for (std::size_t j = i+1; j < minimizers.size() && minimizers[j].itself == minimizers[i].itself; ++j) {++i;}
        }
    }
    assert(none_positions.size() == none_sizes.size());
    std::cout << "lr: " << static_cast<double>(n_maximal) / minimizers.size() * 100 << "%\n";
    std::cout << "l : " << static_cast<double>(left_positions.size()) / minimizers.size() * 100 << "%\n";
    std::cout << "r : " << static_cast<double>(right_or_collision_sizes.size() - colliding_minimizers.size()) / minimizers.size() * 100 << "%\n";
    std::cout << "n : " << static_cast<double>(none_positions.size()) / minimizers.size() * 100 << "%\n";
    std::cout << "ambiguous : " << colliding_minimizers.size() << "/" << n_distinct_minimizers << " (" << static_cast<double>(colliding_minimizers.size()) / n_distinct_minimizers * 100 << ")%\n";

    quartet_wtree wtree;
    wtree.build(wtb);
    /*part 5: Elias-Fano*/
    struct posize_idx_t {
        // Left positions | right + coll sizes | none sizes | none positions
        std::size_t right_coll_sizes_start, none_sizes_start, none_pos_start;
    };
    posize_idx_t index;
    index.right_coll_sizes_start = left_positions.size();
    left_positions.insert(left_positions.end(), right_or_collision_sizes.begin(), right_or_collision_sizes.end());
    right_or_collision_sizes.clear();

    index.none_sizes_start = left_positions.size();
    left_positions.insert(left_positions.end(), none_sizes.begin(), none_sizes.end());
    none_sizes.clear();

    index.none_pos_start = left_positions.size();
    left_positions.insert(left_positions.end(), none_positions.begin(), none_positions.end());
    none_positions.clear();

    pthash::ef_sequence<true> inv_idx;
    inv_idx.encode(left_positions.begin(), left_positions.size());

    /*part 6: build fallback mphf*/
    std::cout << "Part 6: fallback MPHF\n";
    std::sort(colliding_minimizers.begin(), colliding_minimizers.end()); // FIXME sort minimizers for fast search -> find better alternative (hash table)
    std::vector<uint64_t> unbucketable_kmers;
    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) { // reopen input file in order to find ambiguous minimizers and their k-mers
        std::cerr << "Unable to open the input file a second time" << input_filename << "\n";
        return 2;
    }
    seq = kseq_init(fp);
    while(kseq_read(seq) >= 0) {
        // std::cerr << "--------------------------------------------------------------------------" << std::endl;
        std::string contig = std::string(seq->seq.s); //we lose a little bit of efficiency here
        minimizer::get_colliding_kmers<murmurhash2_64>(contig, k, m, seed, false, colliding_minimizers, unbucketable_kmers);
        // std::cerr << "//////////////////////////////////////////////////////////////////////////" << std::endl;
    }
    if (seq) kseq_destroy(seq);
    
    // std::cerr << "colliding k-mers : " << unbucketable_kmers << std::endl;

    pthash_mphf_type kmer_mphf;
    kmer_mphf.build_in_external_memory(unbucketable_kmers.begin(), unbucketable_kmers.size(), mphf_config);

    /*part 7: check for correctness (query + check minimality)*/
    std::cout << "Part 7: check\n";
    pthash::bit_vector_builder population(total_kmers); // bitvector for checking perfection and minimality
    if ((fp = gzopen(input_filename.c_str(), "r")) == NULL) { // reopen input stream once again
        std::cerr << "Unable to open the input file a second time" << input_filename << "\n";
        return 2;
    }
    seq = kseq_init(fp);
    while(kseq_read(seq) >= 0) {
        // std::cerr << "--------------------------------------------------------------------------" << std::endl;
        std::string contig = std::string(seq->seq.s); //we lose a little bit of efficiency here
        for (std::size_t i = 0; i < contig.size() - k + 1; ++i) { // This is NOT streaming unlike the previous walks. FIXME: when making a separate module for query (for the optimized version of the paper)
            uint64_t kmer = debug::string_to_uint64_no_reverse(&contig[i], k);
            debug::triplet_t triplet = debug::compute_minimizer_triplet(kmer, k, m, seed);
            uint64_t mm = triplet.first;
            uint64_t p = triplet.third;
            uint64_t mp_hash = mm_mphf(mm);
            // std::cerr << "mm hash = " << mp_hash << "\n";

            std::pair<MinimizerType, std::size_t> dummy = wtree.rank_of(mp_hash);
            MinimizerType mm_type = dummy.first;
            uint64_t mm_type_rank = dummy.second;
            uint64_t locpres_hash, sk_size;

            switch (mm_type) {
                case LEFT: 
                    locpres_hash = 0; // because in the elias-fano global vector left positions are the left-most block starting at the beginning
                    // std::cerr << "[LEFT] rank = " << mm_type_rank << "\n";
                    locpres_hash += inv_idx.my_access(mm_type_rank); // number of left-KMERS before our bucket
                    // std::cerr << "[LEFT] global shift = " << locpres_hash << ", local rank = " << p << "\n"; 
                    locpres_hash += p; // add local rank
                    break;
                case RIGHT_OR_COLLISION:
                    // std::cerr << "[RIGHT/COLLISION] rank = " << index.right_coll_sizes_start + mm_type_rank << "\n";
                    locpres_hash = inv_idx.my_access(index.right_coll_sizes_start + mm_type_rank); // global shift
                    // std::cerr << "[RIGHT] global shift = " << locpres_hash << ", "; 
                    sk_size = inv_idx.diff(index.right_coll_sizes_start + mm_type_rank);
                    if (sk_size == 0) {
                        // std::cerr << "[COLLISION] rank = " << index.none_pos_start << "\n";
                        locpres_hash = inv_idx.my_access(index.none_pos_start); // prefix sum of all sizes (sizes of collisions are 0)
                        locpres_hash += kmer_mphf(kmer);
                    } else {
                        // std::cerr << "local shift = " << k - m - p << "\n";
                        locpres_hash += k - m - p; // local shift
                    }
                    break;
                case MAXIMAL: // easy case
                    // std::cerr << "[MAXIMAL] rank = " << mm_type_rank << "\n";
                    // std::cerr << "[MAXIMAL] global shift = " << (k-m+1) * mm_type_rank << ", local shift = " << p << "\n";
                    locpres_hash = (k-m+1) * mm_type_rank + p; // all maximal k-mer hashes are < than those of all the other types
                    break;
                case NONE:
                    // std::cerr << "[NONE] rank = " << index.none_sizes_start + mm_type_rank << " = " << index.none_sizes_start << " + " << mm_type_rank << "\n";
                    locpres_hash = inv_idx.my_access(index.none_sizes_start + mm_type_rank); // prefix sum of sizes
                    // std::cerr << "[NONE] global rank = " << locpres_hash << ", ";
                    sk_size = inv_idx.diff(index.none_pos_start + mm_type_rank); // p1 actually
                    locpres_hash += sk_size - p; // position in the first k-mer - actual position = local shift
                    // std::cerr << "[NONE] local rank = " << sk_size - p << ", p1 = " << sk_size << ", p = " << p << "\n";
                    break;
                default:
                    std::cerr << "[Error] Something went wrong with the Wavelet Tree, minimizer type not recognised" << std::endl;
                    return 128;
            }
            if (mm_type != MAXIMAL) locpres_hash += (k-m+1) * n_maximal; // shift of the maximal k-mers
            if (false and mm_type == NONE) {
                std::string explicit_kmer(contig, i, k);
                std::cerr << explicit_kmer << ", ";
                std::cerr << "minimizer = " << mm;
                std::cerr << std::dec << ", type = " << mm_type << ", mm pos = " << p << ", hash = " << locpres_hash << "\n";
            }
            if (locpres_hash > total_kmers) {
                std::cerr << "[Error] overflow : " << locpres_hash << " > " << total_kmers << std::endl;
                return 128;
            } else if (population.get(locpres_hash) == 1) { // Error, we saw a collision
                std::cerr << "[Error] collision at position (hash) : " << locpres_hash << std::endl; 
                return 128;
            } else { // ok
                population.set(locpres_hash);
            }
        }
        // std::cerr << "//////////////////////////////////////////////////////////////////////////" << std::endl;
    }
    if (seq) kseq_destroy(seq);
    bool perfect = true;
    for (std::size_t i = 0; i < total_kmers; ++i) {
        if (!population.get(i)) perfect = false;
    }
    if (!perfect) std::cerr << "[Error] Not all k-mers have been marked by a hash" << std::endl;
    else std::cout << "[Info] Everything is ok\n";
    /*part 8: statistics*/
    std::cerr << "\nPart 8: Statistics ---------------------------------------------------------\n";
    auto mm_mphf_size_bits = mm_mphf.num_bits();
    auto triplet_tree_size_bits = wtree.num_bits();
    auto elias_sequence_size_bits = sizeof(posize_idx_t) * 8 + inv_idx.num_bits();
    auto kmer_mphf_size_bits = kmer_mphf.num_bits();
    auto total_bit_size = mm_mphf_size_bits + triplet_tree_size_bits + elias_sequence_size_bits + kmer_mphf_size_bits;
    std::cout << "Minimizer MPHF size in bits : " << mm_mphf_size_bits << " (" << static_cast<double>(mm_mphf_size_bits) / total_bit_size * 100 << "%)\n";
    std::cout << "\t = " << static_cast<double>(mm_mphf_size_bits) / minimizers.size() << " bits/minimizer\n\n";
    std::cout << "Wavelet tree size in bits : " << triplet_tree_size_bits << " (" << static_cast<double>(triplet_tree_size_bits) / total_bit_size * 100 << "%)\n";
    std::cout << "\t = " << static_cast<double>(triplet_tree_size_bits) / minimizers.size() << " bits/minimizer\n\n";
    std::cout << "Compressed arrays (EF) : " << elias_sequence_size_bits << " (" << static_cast<double>(elias_sequence_size_bits) / total_bit_size * 100 << "%)\n";
    std::cout << "\t = " << static_cast<double>(elias_sequence_size_bits) / left_positions.size() << " bits/offset\n\n";
    std::cout << "Fallback MPHF : " << kmer_mphf_size_bits << " (" << static_cast<double>(kmer_mphf_size_bits) / total_bit_size * 100 << "%)\n";
    std::cout << "\t = " << static_cast<double>(kmer_mphf_size_bits) / unbucketable_kmers.size() << " bits/kmer\n\n";
    std::cout << "Total size in bits : " << total_bit_size << "\n";
    std::cout << "\tequivalent to : " << static_cast<double>(total_bit_size) / total_kmers << " bits/k-mer\n";    
    return 0;
}