#include "../include/ptbb_file_itr.hpp"

extern "C" {
#include "../include/kseq.h"
}

KSEQ_INIT(gzFile, gzread)

namespace lphash {

namespace other {

ptbb_file_itr::ptbb_file_itr() : z(0), k(0), fp(nullptr), sequence_file(nullptr), nbases_since_last_break(0), base_index(0), hn(false), ref_count(nullptr) {};

ptbb_file_itr::ptbb_file_itr(std::string fasta_file, uint64_t kmer_len) : z(0), k(kmer_len), fp(nullptr), nbases_since_last_break(0), base_index(0)
{
    km_shift = 2 * (k - 1);
    km_mask = (static_cast<kmer_t>(1) << (2 * k)) - 1;
    km[0] = 0;
    km[1] = 0;
    if ((fp = gzopen(fasta_file.c_str(), "r")) == NULL) 
        throw std::runtime_error("[ptbb_file_itr] Unable to open file");
    kseq_t* seq = kseq_init(fp);
    sequence_file = seq;
    if (kseq_read(seq) < 0) hn = false;
    else {
        if (seq->seq.l < k) base_index = seq->seq.l;
        hn = true;
        operator++();
    }
    ref_count = (uint64_t*) malloc(sizeof(uint64_t));
    *ref_count = 1;
}

ptbb_file_itr::ptbb_file_itr(const ptbb_file_itr & other)
{
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

void ptbb_file_itr::operator++() 
{
    kseq_t* seq = reinterpret_cast<kseq_t*>(sequence_file);
    while(hn && base_index >= seq->seq.l) {
        if (base_index >= seq->seq.l) { // load next sequence, if any
            base_index = 0;
            nbases_since_last_break = 0;
            auto r = kseq_read(seq);
            if (r < 0) hn = false;
            else if (seq->seq.l < k) base_index = seq->seq.l;
        }
    }
    while(hn && nbases_since_last_break < k && base_index < seq->seq.l) {
        c = constants::seq_nt4_table[static_cast<uint8_t>(seq->seq.s[base_index])];
        if (c < 4) {
            km[0] = (km[0] << 2 | static_cast<kmer_t>(c)) & km_mask;
            km[1] = (km[1] >> 2) | ((static_cast<kmer_t>(3) ^ static_cast<kmer_t>(c)) << km_shift);
            // if (km[0] != km[1]) z = km[0] < km[1] ? 0 : 1;
            ++nbases_since_last_break;
        } else {
            nbases_since_last_break = 0;
        }
        ++base_index;
    }
    if (nbases_since_last_break > 0) --nbases_since_last_break;
}

ptbb_file_itr::~ptbb_file_itr() 
{
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

bool operator==(ptbb_file_itr const& a, ptbb_file_itr const& b)
{
    return a.hn == b.hn;
}
bool operator!=(ptbb_file_itr const& a, ptbb_file_itr const& b)
{
    return !(a == b);
}

}

}