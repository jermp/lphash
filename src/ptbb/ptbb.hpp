#pragma once

#include <zlib.h>

extern "C" {
#include "../../external/kseq.h"
}

KSEQ_INIT(gzFile, gzread)

#include "../../include/constants.hpp"
#include "../../external/BooPHF.hpp"

namespace ptbb {

struct PTHasher {
    typedef pthash::hash128 hash_type;
    static inline pthash::hash128 hash(lphash::kmer_t val, uint64_t seed) {
        return {pthash::MurmurHash2_64(reinterpret_cast<char const*>(&val), sizeof(val), seed),
                pthash::MurmurHash2_64(reinterpret_cast<char const*>(&val), sizeof(val), ~seed)};
    }
};

struct BBHasher {
    uint64_t operator()(lphash::kmer_t val, uint64_t seed = 1234567890) const {
        pthash::hash128 hash{
            pthash::MurmurHash2_64(reinterpret_cast<char const*>(&val), sizeof(val), seed),
            pthash::MurmurHash2_64(reinterpret_cast<char const*>(&val), sizeof(val), ~seed)};
        return hash.mix();
    };
};

typedef pthash::single_phf<PTHasher, pthash::dictionary_dictionary, true> pthash_mphf_t;
typedef boomphf::mphf<lphash::kmer_t, BBHasher> bbhash_mphf_t;

class ptbb_file_itr : std::forward_iterator_tag {
public:
    typedef lphash::kmer_t value_type;
    ptbb_file_itr(const ptbb_file_itr&);
    ptbb_file_itr();
    ptbb_file_itr(std::string fasta_file, uint64_t kmer_len);
    ~ptbb_file_itr();
    inline lphash::kmer_t operator*() const { return km[z]; };
    inline bool has_next() const { return hn; };
    inline void* memory_management() {
        return sequence_file;
    };  // FIXME very bad workaround but ok for now.
    void operator++();

private:
    int c, z;
    uint64_t k;
    gzFile fp;
    void* sequence_file;
    std::size_t nbases_since_last_break;
    std::size_t base_index;
    uint64_t km_shift;
    lphash::kmer_t km_mask;
    std::array<lphash::kmer_t, 2> km;
    bool hn;
    uint64_t* ref_count;
    friend bool operator==(ptbb_file_itr const& a, ptbb_file_itr const& b);
};

bool operator==(ptbb_file_itr const& a, ptbb_file_itr const& b) { return a.hn == b.hn; }
bool operator!=(ptbb_file_itr const& a, ptbb_file_itr const& b) { return !(a == b); }

ptbb_file_itr::ptbb_file_itr()
    : z(0)
    , k(0)
    , fp(nullptr)
    , sequence_file(nullptr)
    , nbases_since_last_break(0)
    , base_index(0)
    , hn(false)
    , ref_count(nullptr){};

ptbb_file_itr::ptbb_file_itr(std::string fasta_file, uint64_t kmer_len)
    : z(0), k(kmer_len), fp(nullptr), nbases_since_last_break(0), base_index(0) {
    km_shift = 2 * (k - 1);
    km_mask = (static_cast<lphash::kmer_t>(1) << (2 * k)) - 1;
    km[0] = 0;
    km[1] = 0;
    if ((fp = gzopen(fasta_file.c_str(), "r")) == NULL)
        throw std::runtime_error("[ptbb_file_itr] Unable to open file");
    kseq_t* seq = kseq_init(fp);
    sequence_file = seq;
    if (kseq_read(seq) < 0)
        hn = false;
    else {
        if (seq->seq.l < k) base_index = seq->seq.l;
        hn = true;
        operator++();
    }
    ref_count = (uint64_t*)malloc(sizeof(uint64_t));
    *ref_count = 1;
}

ptbb_file_itr::ptbb_file_itr(ptbb_file_itr const& other) {
    z = other.z;
    k = other.k;
    fp = other.fp;
    sequence_file = other.sequence_file;
    nbases_since_last_break = other.nbases_since_last_break;
    base_index = other.base_index;
    km_shift = other.km_shift;
    km_mask = other.km_mask;
    km = other.km;
    hn = other.hn;
    ref_count = other.ref_count;
    if (ref_count) ++(*ref_count);
}

void ptbb_file_itr::operator++() {
    kseq_t* seq = reinterpret_cast<kseq_t*>(sequence_file);
    while (hn && base_index >= seq->seq.l) {
        if (base_index >= seq->seq.l) {  // load next sequence, if any
            base_index = 0;
            nbases_since_last_break = 0;
            auto r = kseq_read(seq);
            if (r < 0)
                hn = false;
            else if (seq->seq.l < k)
                base_index = seq->seq.l;
        }
    }
    while (hn && nbases_since_last_break < k && base_index < seq->seq.l) {
        c = lphash::constants::seq_nt4_table[static_cast<uint8_t>(seq->seq.s[base_index])];
        if (c < 4) {
            km[0] = (km[0] << 2 | static_cast<lphash::kmer_t>(c)) & km_mask;
            km[1] = (km[1] >> 2) |
                    ((static_cast<lphash::kmer_t>(3) ^ static_cast<lphash::kmer_t>(c)) << km_shift);
            ++nbases_since_last_break;
        } else {
            nbases_since_last_break = 0;
        }
        ++base_index;
    }
    if (nbases_since_last_break > 0) --nbases_since_last_break;
}

ptbb_file_itr::~ptbb_file_itr() {
    if (ref_count) {
        if ((*ref_count) == 1) {
            if (sequence_file) kseq_destroy(reinterpret_cast<kseq_t*>(sequence_file));
            if (fp) gzclose(fp);
            free(ref_count);
        } else if (*ref_count > 1) {
            --(*ref_count);
        }
    }
}

std::vector<std::string> load_kmers(std::string fasta_file, uint64_t kmer_len) {
    gzFile fp = nullptr;
    kseq_t* seq = nullptr;
    if ((fp = gzopen(fasta_file.c_str(), "r")) == NULL) 
        throw std::runtime_error("Unable to open the input file " + fasta_file + "\n");

    seq = kseq_init(fp);
    std::vector<std::string> kmers;
    while (kseq_read(seq) >= 0) {
        for (std::size_t i = 0; i < seq->seq.l - kmer_len + 1; ++i) {
            kmers.emplace_back(seq->seq.s[i], kmer_len);
        }
    }
    if (seq) kseq_destroy(seq);
    gzclose(fp);
    return kmers;
}

}  // namespace ptbb