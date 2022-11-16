#pragma once

#include <vector>
#include "../include/constants.hpp"
#include "../include/mm_quartet.hpp"
#include "../include/external_memory_vector.hpp"

namespace lphash::minimizer {

template <typename MinimizerHasher>
[[nodiscard]] uint64_t from_string(char const* contig, std::size_t contig_size, uint32_t k,
                                   uint32_t m, uint64_t seed, bool canonical_m_mers,
                                   uint64_t& mm_count,
                                   lphash::external_memory_vector<mm_record_t>& accumulator) {
    std::size_t buf_pos, min_pos;
    mm_quartet_t current;
    uint64_t shift = 2 * (m - 1);
    uint64_t mask = (1ULL << (2 * m)) - 1;
    uint64_t mm[2] = {0, 0};
    uint64_t nbases_since_last_break = 0;
    uint32_t sks = 0, p1 = 0;
    uint64_t kmer_count;
    std::vector<mm_quartet_t> buffer(k - m + 1);
    int c;
    uint8_t z;
    bool find_brand_new_min = false;

    auto update_output = [](decltype(accumulator)& accumulator, mm_quartet_t const& added) {
        // std::cerr << added.itself << " " << unique_mm_id << " " <<
        // static_cast<uint64_t>(added.p1) << " " << static_cast<uint64_t>(added.size) << "\n";
        accumulator.push_back({added.itself, added.id, added.p1, added.size});
    };

    assert(k >= m);
    // std::cerr << "buffer size = " << buffer.size() << "\n";

    buf_pos = 0;
    min_pos = buffer.size();
    kmer_count = 0;
    z = 0;
    for (uint64_t i = 0; i < contig_size; ++i) {
        c = constants::seq_nt4_table[static_cast<uint8_t>(contig[i])];
        current.clear();
        if (c < 4) [[likely]] {
                mm[0] = (mm[0] << 2 | c) & mask;            /* forward m-mer */
                mm[1] = (mm[1] >> 2) | (3ULL ^ c) << shift; /* reverse m-mer */
                if (canonical_m_mers && mm[0] != mm[1])
                    z = mm[0] < mm[1] ? 0
                                      : 1;  // strand, if symmetric k-mer then use previous strand
                ++nbases_since_last_break;
                if (nbases_since_last_break >= m) {
                    current.itself = mm[z];
                    // std::cerr << current.itself << "\n";
                    current.hash = MinimizerHasher::hash(mm[z], seed)
                                       .first();  // insert new hash inside buffer
                    current.p1 =
                        i - m + 1;  // FIXME this is NOT the position inside the super-k-mer!
                    current.id = mm_count++;
                    if (nbases_since_last_break == k) ++kmer_count;
                    if (nbases_since_last_break == k + 1)
                        [[unlikely]] {  // have seen the first window after a break, time to search
                                        // for the minimum note that the current m-mer is checked by
                                        // the next if
                            min_pos = p1 = 0;
                            for (std::size_t j = 0; j < buffer.size(); ++j) {
                                if (buffer[j].hash < buffer[min_pos].hash) {
                                    min_pos = j;
                                    p1 = min_pos;
                                    // std::cerr << "[new min] min_pos = " << min_pos << " " <<
                                    // buffer[min_pos].itself << "\n";
                                }
                            }
                            // std::cerr << buffer << std::endl;
                            // std::cerr << "min_pos = " << min_pos << " " << buffer[min_pos].itself
                            // << "\n";
                            sks = 1;  // min_pos + 1; // number of k-mers after a break is 1
                        }
                    // std::cerr << "partial super-k-mer length: " << sks << "\n";
                    if (nbases_since_last_break >= k + 1)
                        [[likely]] {  // time to update the minimum, if necessary
                            // std::cerr << buf_pos << "\n";
                            assert(sks != 0);
                            assert(sks <= k - m + 1);
                            if (((buf_pos) % buffer.size()) ==
                                min_pos) {  // old minimum outside window
                                buffer[min_pos].p1 = p1;
                                buffer[min_pos].size = sks;
                                update_output(accumulator,
                                              buffer[min_pos]);  // we save the old minimum, length
                                                                 // on the right is k by definition
                                sks = 0;
                                find_brand_new_min = true;  // also update p1
                            } else if (current.hash < buffer[min_pos].hash) {
                                buffer[min_pos].p1 = p1;
                                buffer[min_pos].size = sks;
                                update_output(accumulator, buffer[min_pos]);  // new minimum
                                sks = 0;
                                p1 = k - m;
                                min_pos = buf_pos;  // actual update is outside if
                            }
                            ++sks;
                            ++kmer_count;
                        }
                    buffer[buf_pos++] = current;
                    buf_pos %= buffer.size();  // circular buffer
                    if (find_brand_new_min) {  // find new minimum if the old one dropped out the
                                               // window
                        find_brand_new_min = false;
                        min_pos = buf_pos;
                        p1 = 0;
                        uint32_t tmp = 1;
                        for (std::size_t j = (buf_pos + 1) % buffer.size(); j < buffer.size();
                             ++j) {
                            // std::cerr << "2-" << buffer[min_pos].itself << " " <<
                            // buffer[j].itself << " " << tmp << "\n";
                            if (buffer[min_pos].hash > buffer[j].hash) {
                                min_pos = j;
                                p1 = tmp;
                            }
                            ++tmp;
                        }
                        // std::cerr << "----- " << tmp << std::endl;
                        for (std::size_t j = 0; j <= buf_pos; ++j) {
                            // std::cerr << "2-" << buffer[min_pos].itself << " " <<
                            // buffer[j].itself << " " << tmp << "\n";
                            if (buffer[min_pos].hash > buffer[j].hash) {
                                min_pos = j;
                                p1 = tmp;
                            }
                            ++tmp;
                        }
                        // std::cerr << "[after outside] buffer[" << min_pos << "] = " <<
                        // buffer[min_pos] << "\n";
                    }
                }
            }
        else
            [[unlikely]] {
                nbases_since_last_break = 0;
                if (min_pos < buffer.size()) {
                    buffer[min_pos].p1 = p1;
                    buffer[min_pos].size = sks;
                    update_output(accumulator,
                                  buffer[min_pos]);  // push current minimum if available
                }
                sks = 0;  // impossible value, wait for reinitialization of the first window
                min_pos = buffer.size();
                buf_pos = 0;  // we always restart at the beginning of the buffer -> this allows to
                              // use min_pos as the position of the minimizer inside the first k-mer
            }
    }
    if (nbases_since_last_break == k) {  // contig.length == 1
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
        update_output(accumulator, buffer[min_pos]);  // push last minimum if available
        sks = 1;
    }
    return kmer_count;
}

template <typename MinimizerHasher>
void get_colliding_kmers(char const* contig, std::size_t contig_size, uint32_t k, uint32_t m,
                         uint64_t seed, bool canonical_m_mers,
                         external_memory_vector<uint64_t>::const_iterator& itr,
                         external_memory_vector<uint64_t>::const_iterator& stop, uint64_t& mm_count,
                         external_memory_vector<kmer_t, false>& accumulator) {
    std::vector<mm_record_t> mm_buffer(k - m + 1);
    std::vector<kmer_t> km_buffer;
    std::size_t mm_buf_pos = 0, min_pos = mm_buffer.size();
    mm_record_t current;
    uint64_t mm_shift = 2 * (m - 1);
    uint64_t mm_mask = (1ULL << (2 * m)) - 1;
    uint64_t km_shift = 2 * (k - 1);
    kmer_t km_mask = (static_cast<kmer_t>(1) << (2 * k)) - 1;
    uint64_t mm[2] = {0, 0};
    kmer_t km[2] = {0, 0};
    uint64_t nbases_since_last_break = 0;
    uint32_t sks = 0;
    uint8_t z = 0;
    bool find_brand_new_min = false;
    int c;
    assert(k >= m);

    auto update_output = [&accumulator](std::vector<kmer_t> const& toadd) {
        for (auto kmer : toadd) accumulator.push_back(kmer);
        // if (statistics.count(toadd.size()) == 0) {
        //     statistics[toadd.size()] = 1ULL;
        // } else {
        //     ++statistics[toadd.size()];
        // }
    };

    km_buffer.reserve(2 * k - m);
    for (uint64_t i = 0; i < contig_size; ++i) {
        c = constants::seq_nt4_table[static_cast<uint8_t>(contig[i])];
        if (c < 4) [[likely]] {
                mm[0] = (mm[0] << 2 | c) & mm_mask;            /* forward k-mer */
                mm[1] = (mm[1] >> 2) | (3ULL ^ c) << mm_shift; /* reverse k-mer */
                km[0] = (km[0] << 2 | static_cast<kmer_t>(c)) & km_mask;
                km[1] =
                    (km[1] >> 2) | ((static_cast<kmer_t>(3) ^ static_cast<kmer_t>(c)) << km_shift);
                if (canonical_m_mers && mm[0] != mm[1])
                    z = mm[0] < mm[1] ? 0
                                      : 1;  // strand, if symmetric k-mer then use previous strand
                ++nbases_since_last_break;

                if (nbases_since_last_break >= m) {
                    // current.itself = mm[z];
                    current.itself = MinimizerHasher::hash(mm[z], seed)
                                         .first();  // here we use itself as hash field since we are
                                                    // not interested in minimizers
                    current.id = mm_count++;
                    // std::cerr << current.id << " " << *itr << "\n";
                    // assert(current.id <= *itr + 1);
                    if (nbases_since_last_break == k + 1)
                        [[unlikely]] {  // we have seen the first window after a break, time to
                                        // search for the minimum
                            min_pos = 0;
                            for (std::size_t j = 0; j < mm_buffer.size(); ++j) {
                                if (mm_buffer[j].itself < mm_buffer[min_pos].itself) min_pos = j;
                            }
                            sks = 1;  // number of k-mers after a break is 1
                        }
                    if (nbases_since_last_break >= k + 1)
                        [[likely]] {  // time to update the minimum, if necessary
                            // std::cerr << "current minimizer id = " << mm_buffer[min_pos].id <<
                            // "\n";
                            assert(sks != 0);
                            assert(sks <= k - m + 1);
                            if (((mm_buf_pos) % mm_buffer.size()) == min_pos ||
                                current.itself < mm_buffer[min_pos].itself) {  // update min
                                assert(sks == km_buffer.size());
                                if (itr != stop && *itr == mm_buffer[min_pos].id) {
                                    assert(*itr >= mm_buffer[min_pos].id);
                                    update_output(
                                        km_buffer);  // we save all k-mers in the super-k-mer
                                    ++itr;
                                }
                                km_buffer.clear();
                                if (((mm_buf_pos) % mm_buffer.size()) == min_pos) {
                                    find_brand_new_min = true;  // old minimum outside window
                                } else if (current.itself < mm_buffer[min_pos].itself) {
                                    min_pos =
                                        mm_buf_pos;  // new minimum, actual update is outside if
                                }
                                sks = 0;
                            }
                            ++sks;
                        }

                    mm_buffer[mm_buf_pos++] = current;
                    mm_buf_pos %= mm_buffer.size();  // circular buffer
                    if (nbases_since_last_break >= k)
                        km_buffer.push_back(km[z]);  // put k-mer into current super-k-mer
                    if (find_brand_new_min) {  // find new minimum if the old one dropped out the
                                               // window
                        find_brand_new_min = false;
                        min_pos = mm_buf_pos;
                        for (std::size_t j = (mm_buf_pos + 1) % mm_buffer.size();
                             j < mm_buffer.size(); ++j) {
                            if (mm_buffer[min_pos].itself > mm_buffer[j].itself) min_pos = j;
                        }
                        for (std::size_t j = 0; j <= mm_buf_pos; ++j) {
                            if (mm_buffer[min_pos].itself > mm_buffer[j].itself) min_pos = j;
                        }
                    }
                }
            }
        else
            [[unlikely]] {
                nbases_since_last_break = 0;
                if (min_pos <
                    mm_buffer
                        .size())  // conditions && itr != stop && *itr == mm_count are superfluous
                {
                    assert(sks == km_buffer.size());
                    if (itr != stop && *itr == mm_buffer[min_pos].id) {
                        assert(*itr >= mm_buffer[min_pos].id);
                        update_output(km_buffer);
                        ++itr;
                    }
                }
                km_buffer.clear();
                min_pos = mm_buffer.size();
                sks = 0;  // impossible value, wait for reinitialization of the first window
                mm_buf_pos =
                    0;  // we always restart at the beginning of the buffer -> this allows to use
                        // min_pos as the position of the minimizer inside the first k-mer
            }
    }
    if (nbases_since_last_break == k) {  // contig.length == 1
        min_pos = 0;
        sks = 1;
        for (std::size_t j = 0; j < mm_buffer.size(); ++j) {
            if (mm_buffer[j].itself < mm_buffer[min_pos].itself) { min_pos = j; }
        }
    }
    if (min_pos <
        mm_buffer.size())  // conditions && itr != stop && *itr == mm_count are superfluous
    {
        assert(sks == km_buffer.size());
        if (itr != stop && *itr == mm_buffer[min_pos].id) {
            assert(*itr >= mm_buffer[min_pos].id);
            update_output(km_buffer);
            ++itr;
        }
    }
}

std::pair<external_memory_vector<mm_triplet_t, false>, external_memory_vector<uint64_t>> classify(
    external_memory_vector<mm_record_t>& minimizers, uint8_t max_memory, std::string tmp_dirname);

std::pair<external_memory_vector<mm_triplet_t, false>, external_memory_vector<uint64_t>> classify(
    external_memory_vector<mm_record_t>&& minimizers, uint8_t max_memory, std::string tmp_dirname);

}  // namespace lphash::minimizer