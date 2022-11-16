#include "../include/external_memory_vector.hpp"
#include "../include/constants.hpp"

using namespace lphash;

int main() {
    std::cout << "sizeof(record) = " << sizeof(mm_record_t) << "\n";
    external_memory_vector<uint64_t, true> vec64sorted(100, ".", "test_sev");
    for (std::size_t i = 10; i > 0; --i) vec64sorted.push_back(i);
    for (auto itr = vec64sorted.cbegin(); itr != vec64sorted.cend(); ++itr) std::cout << *itr << "\n";

    external_memory_vector<uint64_t, false> vec64unsorted(100, ".", "test_sev");
    for (std::size_t i = 10; i > 0; --i) vec64unsorted.push_back(i);
    for (auto itr = vec64unsorted.cbegin(); itr != vec64unsorted.cend(); ++itr) std::cout << *itr << "\n";
    // external_memory_vector<mm_record_t> vecmmr(
    //     1000, [](mm_record_t const& a, mm_record_t const& b) { return a.itself < b.itself; }, ".",
    //     "test_sev");

    // mm_record_t record;
    // for (std::size_t i = 0; i < 10; ++i) {
    //     record.itself = i;
    //     record.id = i;
    //     record.p1 = 0;
    //     record.size = 0;
    //     vecmmr.push_back(record);
    // }
    // for (auto itr = vecmmr.cbegin(); itr != vecmmr.cend(); ++itr) {
    //     record = *itr;
    //     std::cout << record.itself << "\n";
    // }
    return 0;
}