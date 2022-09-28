#include <array>
#include <vector>
#include <ostream>
#include "../include/constants.hpp"

#include "../include/prettyprint.hpp"

namespace lphash {

struct mmp_t {
    mmp_t() : hash(0), p1(0), size(0) {};
    void clear() noexcept 
    {
        itself = 0; // AAA...A
        hash = std::numeric_limits<decltype(hash)>::max(); 
        p1 = std::numeric_limits<decltype(p1)>::max();
        size = std::numeric_limits<decltype(size)>::max();
    };

    friend std::ostream& operator<<(std::ostream& os, mmp_t const & other) {
        return os << other.itself << " " << other.hash << ' ' << other.p1 << ' ' << other.size;
    }

    uint64_t hash; // minimizer hash
    uint64_t itself; // 2-bit minimizer itself
    uint32_t p1; // position inside first k-mer of the super-k-mer
    uint32_t size; // size (number of k-mers) in the super-k-mer
};

namespace minimizer {

template <typename MinimizerHasher>
[[nodiscard]] 
uint64_t from_string(std::string const & contig, uint32_t k, uint32_t m, uint64_t seed, bool canonical_m_mers, std::vector<mmp_t> & accumulator)
{
    std::size_t buf_pos, min_pos;
    mmp_t current; 
    uint64_t shift = 2 * (m - 1);
    uint64_t mask = (1ULL << (2 * m)) - 1;
    uint64_t mm[2] = {0, 0};
    uint64_t nbases_since_last_break = 0;
    uint32_t sks = 0, p1 = 0;
    uint64_t kmer_count;
    std::vector<mmp_t> buffer(k-m+1);
    int c;
    uint8_t z;
    bool find_brand_new_min = false;

    auto update_output = [](std::vector<mmp_t> & accumulator, mmp_t const & added, std::string const & msg) {
        // std::cerr << "[" << msg << "] to_be_added = {" << added << "}" << std::endl;
        accumulator.push_back(added);
    };

    assert(k >= m);

    buf_pos = 0;
    min_pos = buffer.size();
    // std::cerr << "buffer size = " << buffer.size() << "\n";
    kmer_count = 0;
    z = 0;
    for (uint64_t i = 0; i < contig.size(); ++i) 
    {
        c = constants::seq_nt4_table[static_cast<uint8_t>(contig[i])];
        current.clear();
        if (c < 4) [[likely]] {
			mm[0] = (mm[0] << 2 | c) & mask;          /* forward m-mer */
			mm[1] = (mm[1] >> 2) | (3ULL^c) << shift; /* reverse m-mer */
			// if (canonical_m_mers and mm[0] != mm[1]) z = mm[0] < mm[1] ? 0 : 1; // strand, if symmetric k-mer then use previous strand
			++nbases_since_last_break;
			if (nbases_since_last_break >= m) {
                current.itself = mm[z];
                // std::cerr << current.itself << "\n";
                current.hash = MinimizerHasher::hash(mm[z], seed);// insert new hash inside buffer
                current.p1 = i-m+1; // FIXME this is NOT the position inside the super-k-mer!
                if (nbases_since_last_break == k) ++kmer_count;
                if (nbases_since_last_break == k + 1) [[unlikely]] { // have seen the first window after a break, time to search for the minimum
                    // note that the current m-mer is checked by the next if
                    min_pos = p1 = 0;
                    for (std::size_t j = 0; j < buffer.size(); ++j) {
                        if (buffer[j].hash < buffer[min_pos].hash) {
                            min_pos = j;
                            p1 = min_pos;
                            // std::cerr << "[new min] min_pos = " << min_pos << " " << buffer[min_pos].itself << "\n";
                        }
                    }
                    //std::cerr << buffer << std::endl;
                    // std::cerr << "min_pos = " << min_pos << " " << buffer[min_pos].itself << "\n";
                    sks = 1;//min_pos + 1; // number of k-mers after a break is 1
                    // ++kmer_count;
                }
                // std::cerr << "partial super-k-mer length: " << sks << "\n";
                if (nbases_since_last_break >= k + 1) [[likely]] { // time to update the minimum, if necessary
                    // std::cerr << buf_pos << "\n";
                    assert(sks != 0);
                    assert(sks <= k-m+1);
                    if (((buf_pos) % buffer.size()) == min_pos) { // old minimum outside window
                        buffer[min_pos].p1 = p1;
                        buffer[min_pos].size = sks;
                        update_output(accumulator, buffer[min_pos], "outside"); // we save the old minimum, length on the right is k by definition
                        sks = 0;
                        find_brand_new_min = true; // also update p1
                    } else if (current.hash < buffer[min_pos].hash) {
                        // uint32_t old_idx = buffer[min_pos].p1;
                        buffer[min_pos].p1 = p1;
                        buffer[min_pos].size = sks;
                        update_output(accumulator, buffer[min_pos], "new"); // new minimum
                        sks = 0;
                        p1 = k-m;
                        min_pos = buf_pos; // actual update is outside if
                    }
                    ++sks;
                    ++kmer_count;
                }
                buffer[buf_pos++] = current;
                buf_pos %= buffer.size(); // circular buffer
                if (find_brand_new_min) { // find new minimum if the old one dropped out the window
                    find_brand_new_min = false;
                    min_pos = buf_pos;
                    p1 = 0;
                    uint32_t tmp = 1;
                    for (std::size_t j = (buf_pos + 1) % buffer.size(); j < buffer.size(); ++j) {
                        // std::cerr << "2-" << buffer[min_pos].itself << " " << buffer[j].itself << " " << tmp << "\n";
                        if (buffer[min_pos].hash > buffer[j].hash) {
                            min_pos = j;
                            p1 = tmp;
                            
                        }
                        ++tmp;
                    }
                    // std::cerr << "----- " << tmp << std::endl;
                    for (std::size_t j = 0; j <= buf_pos; ++j) {
                        // std::cerr << "2-" << buffer[min_pos].itself << " " << buffer[j].itself << " " << tmp << "\n";
                        if (buffer[min_pos].hash > buffer[j].hash) {
                            min_pos = j;
                            p1 = tmp;
                            
                        }
                        ++tmp;
                    }
                    // std::cerr << "[after outside] buffer[" << min_pos << "] = " << buffer[min_pos] << "\n";
                }
            }
        } else [[unlikely]] {
            nbases_since_last_break = 0;
            if (min_pos < buffer.size()) {
                buffer[min_pos].p1 = p1;
                buffer[min_pos].size = sks;
                update_output(accumulator, buffer[min_pos], "last"); //push current minimum if available
            }
            sks = 0; // impossible value, wait for reinitialization of the first window
            min_pos = buffer.size();
            buf_pos = 0; // we always restart at the beginning of the buffer -> this allows to use min_pos as the position of the minimizer inside the first k-mer
        }   
    }
    if (nbases_since_last_break == k) { // contig.length == 1
        min_pos = p1 = 0;
        sks = 1;
        for (std::size_t j = 0; j < buffer.size(); ++j) {
            if (buffer[j].hash < buffer[min_pos].hash) {
                min_pos = j;
                p1 = min_pos;
            }
        }
    }
    if (min_pos < buffer.size()) {
        buffer[min_pos].p1 = p1;
        buffer[min_pos].size = sks;
        update_output(accumulator, buffer[min_pos], "last"); //push last minimum if available  
        sks = 1;
    }
    return kmer_count;
}

template <typename MinimizerHasher, typename KMerType = uint64_t>
void get_colliding_kmers(std::string const & contig, uint32_t k, uint32_t m, uint64_t seed, bool canonical_m_mers, std::vector<uint64_t> const& colliding_minimizers, std::vector<KMerType> & accumulator)
{
    typedef std::pair<uint64_t, uint64_t> mm_pair_t;
    std::vector<mm_pair_t> mm_buffer(k-m+1);
    std::vector<KMerType> km_buffer;
    std::size_t mm_buf_pos = 0, min_pos = mm_buffer.size();
    mm_pair_t current;
    uint64_t mm_shift = 2 * (m - 1);
    uint64_t mm_mask = (1ULL << (2 * m)) - 1;
    uint64_t km_shift = 2 * (k - 1);
    uint64_t km_mask = (1ULL << (2 * k)) - 1;
    uint64_t mm[2] = {0, 0};
    KMerType km[2] = {0, 0};
    uint64_t nbases_since_last_break = 0;
    uint32_t sks = 0;
    uint8_t z = 0;
    bool find_brand_new_min = false;
    int c;
    assert(k >= m);

    essentials::timer_type timer;

    auto update_output = [&](std::vector<KMerType> const& toadd, std::vector<KMerType> & accumulator) {
        accumulator.insert(accumulator.end(), toadd.begin(), toadd.end());
    };
    km_buffer.reserve(2*k-m);
    for (uint64_t i = 0; i < contig.size(); ++i) 
    {
        c = constants::seq_nt4_table[static_cast<uint8_t>(contig[i])];
        if (c < 4) [[likely]] {
			mm[0] = (mm[0] << 2 | c) & mm_mask;          /* forward k-mer */
			mm[1] = (mm[1] >> 2) | (3ULL^c) << mm_shift; /* reverse k-mer */
            km[0] = (km[0] << 2 | c) & km_mask;
            km[1] = (km[1] >> 2) | (3ULL^c) << km_shift;
			// if (canonical_m_mers && mm[0] != mm[1]) z = mm[0] < mm[1] ? 0 : 1; // strand, if symmetric k-mer then use previous strand
			++nbases_since_last_break;

			if (nbases_since_last_break >= m) {
                current.first = mm[z];
                current.second = MinimizerHasher::hash(mm[z], seed);// insert new hash inside buffer
                if (nbases_since_last_break == k + 1) [[unlikely]] { // we have seen the first window after a break, time to search for the minimum
                    min_pos = 0;
                    for (std::size_t j = 0; j < mm_buffer.size(); ++j) {
                        if (mm_buffer[j].second < mm_buffer[min_pos].second) 
                            min_pos = j;
                    }
                    sks = 1; // number of k-mers after a break is 1
                }
                // std::cerr << "partial super-k-mer length: " << sks << "\n";
                if (nbases_since_last_break >= k + 1) [[likely]] { // time to update the minimum, if necessary
                    // std::cerr << buf_pos << "\n";
                    assert(sks != 0);
                    assert(sks <= k-m+1);
                    // std::cerr << "super-k-mer length = " << sks << ", super-k-mer window length = " << km_buffer.size() << std::endl;
                    if (((mm_buf_pos) % mm_buffer.size()) == min_pos || current.second < mm_buffer[min_pos].second) { // update min
                        if (std::binary_search(colliding_minimizers.begin(), colliding_minimizers.end(), mm_buffer[min_pos].first)) {
                            // std::cerr << "[update] super-k-mer length = " << sks << ", super-k-mer window length = " << km_buffer.size() << std::endl;
                            assert(sks == km_buffer.size());
                            update_output(km_buffer, accumulator); // we save all k-mers in the super-k-mer
                        }
                        km_buffer.clear();
                        if (((mm_buf_pos) % mm_buffer.size()) == min_pos) find_brand_new_min = true; // old minimum outside window
                        else if (current.second < mm_buffer[min_pos].second) min_pos = mm_buf_pos; // new minimum, actual update is outside if
                        sks = 0;
                    }
                    ++sks;
                }

                mm_buffer[mm_buf_pos++] = current; 
                mm_buf_pos %= mm_buffer.size(); // circular buffer
                if (nbases_since_last_break >= k) km_buffer.push_back(km[z]); // put k-mer into current super-k-mer

                if (find_brand_new_min) { // find new minimum if the old one dropped out the window
                    find_brand_new_min = false;
                    min_pos = mm_buf_pos;
                    for (std::size_t j = (mm_buf_pos + 1) % mm_buffer.size(); j < mm_buffer.size(); ++j) 
                        if (mm_buffer[min_pos].second > mm_buffer[j].second) min_pos = j;
                    for (std::size_t j = 0; j <= mm_buf_pos; ++j) 
                        if (mm_buffer[min_pos].second > mm_buffer[j].second) min_pos = j;
                }
            }
        } else [[unlikely]] {
            nbases_since_last_break = 0;
            if (min_pos < mm_buffer.size() && std::binary_search(colliding_minimizers.begin(), colliding_minimizers.end(), mm_buffer[min_pos].first)) {
                // std::cerr << "[last after N] super-k-mer length = " << sks << ", super-k-mer window length = " << km_buffer.size() << std::endl;
                assert(sks == km_buffer.size());
                update_output(km_buffer, accumulator); // we save all k-mers in the super-k-mer
            }
            km_buffer.clear();
            min_pos = mm_buffer.size();
            sks = 0; // impossible value, wait for reinitialization of the first window
            mm_buf_pos = 0; // we always restart at the beginning of the buffer -> this allows to use min_pos as the position of the minimizer inside the first k-mer
        }   
    }
    if (nbases_since_last_break == k) { // contig.length == 1
        min_pos = 0;
        sks = 1;
        for (std::size_t j = 0; j < mm_buffer.size(); ++j) {
            if (mm_buffer[j].second < mm_buffer[min_pos].second) {
                min_pos = j;
            }
        }
    }
    if (min_pos < mm_buffer.size() && std::binary_search(colliding_minimizers.begin(), colliding_minimizers.end(), mm_buffer[min_pos].first)) {
        // std::cerr << "[very last] super-k-mer length = " << sks << ", super-k-mer window length = " << km_buffer.size() << std::endl;
        assert(sks == km_buffer.size());
        update_output(km_buffer, accumulator); // we save all k-mers in the super-k-mer
    }
    // std::cerr << "buffer size = " << km_buffer.size() << "\n";
}

} // namespace minimizers
} // namespace lphash