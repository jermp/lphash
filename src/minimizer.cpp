#include "../include/minimizer.hpp"

namespace lphash::minimizer {

std::pair<external_memory_vector<mm_triplet_t>, external_memory_vector<uint64_t>> classify(
    external_memory_vector<mm_record_t>& minimizers, uint8_t max_memory, std::string tmp_dirname) {
    auto start = minimizers.cbegin();  // this forces the remaining buffer to be written to disk
    auto stop = minimizers.cend();
    uint64_t colliding_mm_size_estimate =
        static_cast<uint64_t>(static_cast<double>(minimizers.size()) * 0.01 * sizeof(uint64_t));
    colliding_mm_size_estimate =
        colliding_mm_size_estimate < 4000000 ? 4000000 : colliding_mm_size_estimate;
    uint64_t unique_minimizer_mm_size_estimate =
        uint64_t(max_memory) * essentials::GB - colliding_mm_size_estimate;
    external_memory_vector<mm_triplet_t> unique_minimizers(
        unique_minimizer_mm_size_estimate,
        []([[maybe_unused]] mm_triplet_t const& a, [[maybe_unused]] mm_triplet_t const& b) {
            return false;
        },
        tmp_dirname, get_group_id());
    external_memory_vector<uint64_t> colliding_minimizer_ids(
        colliding_mm_size_estimate, [](uint64_t a, uint64_t b) { return a < b; }, tmp_dirname,
        get_group_id());

    mm_record_t prev;
    prev.size = prev.p1 = 0;
    while (start != stop) {
        if (prev.size != 0) {
            if (prev.itself == (*start).itself) {
                prev.p1 = prev.size = 0;
                unique_minimizers.push_back({prev.itself, prev.p1, prev.size});
                colliding_minimizer_ids.push_back(prev.id);
                while (start != stop && (*start).itself == prev.itself) {
                    colliding_minimizer_ids.push_back((*start).id);
                    ++start;
                }
            } else {
                unique_minimizers.push_back({prev.itself, prev.p1, prev.size});
                prev = *start;
                ++start;
            }
        } else {
            prev = *start;
            ++start;
        }
    }
    if (prev.size) unique_minimizers.push_back({prev.itself, prev.p1, prev.size});

    return std::make_pair(std::move(unique_minimizers), std::move(colliding_minimizer_ids));
}

std::pair<external_memory_vector<mm_triplet_t>, external_memory_vector<uint64_t>> classify(
    external_memory_vector<mm_record_t>&& minimizers, uint8_t max_memory, std::string tmp_dirname) {
    return classify(minimizers, max_memory, tmp_dirname);
}

}  // namespace lphash::minimizer