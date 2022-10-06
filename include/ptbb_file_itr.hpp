#ifndef PTHASH_FILE_ITR_HPP
#define PTHASH_FILE_ITR_HPP

#include <zlib.h>
#include "constants.hpp"

namespace lphash {

namespace other {

class ptbb_file_itr : std::forward_iterator_tag 
{
    public:
        typedef kmer_t value_type;
        ptbb_file_itr(const ptbb_file_itr&);
        ptbb_file_itr();
        ptbb_file_itr(std::string fasta_file, uint64_t kmer_len);
        ~ptbb_file_itr();
        inline kmer_t operator*() const { return km[z]; };
        inline bool has_next() const { return hn; };
        inline void* memory_management() {return sequence_file;}; // FIXME very bad workaround but ok for now.
        void operator++();

    private:
        int c, z;
        uint64_t k;
        gzFile fp;
        void* sequence_file;
        std::size_t nbases_since_last_break;
        std::size_t base_index;
        uint64_t km_shift;
        kmer_t km_mask;
        std::array<kmer_t, 2> km;
        bool hn;
        uint64_t* ref_count;
        //std::string filename;
        friend bool operator==(ptbb_file_itr const& a, ptbb_file_itr const& b);
};

bool operator==(ptbb_file_itr const& a, ptbb_file_itr const& b); // very ugly work-around for bbhash
bool operator!=(ptbb_file_itr const& a, ptbb_file_itr const& b);

template <typename KeyType>
struct BBHasher {
    uint64_t operator() (const KeyType val, uint64_t seed = 1234567890) const {
        uint64_t hval = hash64::hash(val, seed).first();
        // std::cerr << "hval = " << hval << "\n";
        return hval;
    };
};

}

}

#endif // PTHASH_FILE_ITR_HPP