#ifndef SORTED_EXTERNAL_VECTOR_HPP
#define SORTED_EXTERNAL_VECTOR_HPP

#include <string>
#include <sstream>
#include <functional>
#include "constants.hpp"

namespace lphash {

template <typename ItrType>
class append_iterator {
		public:
	typedef typename ItrType::value_type value_type;
	append_iterator(std::vector<std::pair<ItrType, ItrType>>& start_end) : i(0), iterators(start_end) {}

	void operator++() {
		++(iterators[i].first);
		if (iterators[i].first == iterators[i].second) ++i;
	}

	typename ItrType::value_type operator*() const { return *(iterators[i].first); }

	bool has_next() const { // FIXME This is not very C++-ish, rewrite the class as a dummy container with methods begin() and end()
		if (i == iterators.size()) return false;
		return true;
	}

		private:
	uint64_t i;
	std::vector<std::pair<ItrType, ItrType>>& iterators;
};

template <typename T>
class sorted_external_vector {
		public:
	class mm_loader_iterator {
			public:
		mm_loader_iterator(uint8_t const* begin, uint8_t const* end) : m_begin(begin), m_end(end) {}
		void operator++() { m_begin += sizeof(T); }
		bool has_next() const { return m_begin != m_end; }
		T const& operator*() const { return *reinterpret_cast<T const*>(m_begin); }

			private:
		uint8_t const* m_begin;
		uint8_t const* m_end;
	};

	class const_iterator : std::forward_iterator_tag {
			public:
		typedef T value_type;
		const_iterator();
		const_iterator(sorted_external_vector<T> const* vec);
		T const& operator*() const;
		void operator++();
		bool operator==(const_iterator const& other) const;
		bool operator!=(const_iterator const& other) const;

			private:
		sorted_external_vector<T> const* v;
		std::vector<mm_loader_iterator> m_iterators;
		std::vector<uint32_t> m_idx_heap;
		std::vector<mm::file_source<uint8_t>> m_mm_files;
		std::function<bool(uint32_t, uint32_t)> heap_idx_comparator;
		void advance_heap_head();
	};

        sorted_external_vector(uint64_t available_space_bytes, std::function<bool(T const&, T const&)> cmp, std::string tmp_dir, std::string name = "");
        sorted_external_vector(uint64_t available_space_bytes, std::string tmp_dir, std::string name = "");
        sorted_external_vector(sorted_external_vector &&) = default;
        void push_back(T const& elem);
        const_iterator cbegin();
        const_iterator cend() const;
        std::size_t size() const;
        ~sorted_external_vector();

		private:
	void init(uint64_t available_space_bytes);
	std::function<bool(T const&, T const&)> m_sorter;
	std::size_t m_buffer_size;
	std::size_t m_total_elems;
	std::string m_tmp_dirname;
	std::string m_prefix;
	std::vector<std::string> m_tmp_files;
	std::vector<T> m_buffer;

	void sort_and_flush();
	std::string get_tmp_output_filename(uint64_t id) const;
};

template <typename T>
sorted_external_vector<T>::sorted_external_vector(uint64_t available_space_bytes, std::function<bool(T const&, T const&)> cmp, std::string tmp_dir, std::string name)
	: m_sorter(cmp), m_total_elems(0), m_tmp_dirname(tmp_dir), m_prefix(name) {
	init(available_space_bytes);
}

template <typename T>
sorted_external_vector<T>::sorted_external_vector(uint64_t available_space_bytes, std::string tmp_dir, std::string name)
	: m_sorter([](T const& a, T const& b) { return a < b; }), m_total_elems(0), m_tmp_dirname(tmp_dir), m_prefix(name) {
	init(available_space_bytes);
}

template <typename T>
void sorted_external_vector<T>::init(uint64_t available_space_bytes) {
	m_buffer_size = available_space_bytes / sizeof(T) + 1;
	m_buffer.reserve(m_buffer_size);
}

template <typename T>
void sorted_external_vector<T>::push_back(T const& elem) {
	// for optimal memory management in the general case one should try to reload the last tmp file if it isn't full
	m_buffer.reserve(m_buffer_size);
	m_buffer.push_back(elem);
	++m_total_elems;
	// std::cerr << m_buffer.size() << " >= " << m_buffer_size << std::endl;
	if (m_buffer.size() >= m_buffer_size) {
		// std::cerr << "buffer full: " << m_buffer.size() << std::endl;
		sort_and_flush();
	}
}

template <typename T>
typename sorted_external_vector<T>::const_iterator sorted_external_vector<T>::cbegin()
{
    if (m_buffer.size() != 0) {
        sort_and_flush();
    }
    m_buffer.clear();
    m_buffer.shrink_to_fit();
    return const_iterator(this);
}

template <typename T>
typename sorted_external_vector<T>::const_iterator sorted_external_vector<T>::cend() const {
	return const_iterator();
}

template <typename T>
std::size_t sorted_external_vector<T>::size() const {
	return m_total_elems;
}

template <typename T>
sorted_external_vector<T>::~sorted_external_vector() {
	for (auto tmp : m_tmp_files) std::remove(tmp.c_str());
}

template <typename T>
void sorted_external_vector<T>::sort_and_flush() {
	std::sort(m_buffer.begin(), m_buffer.end(), m_sorter);
	m_tmp_files.push_back(get_tmp_output_filename(m_tmp_files.size()));
	std::ofstream out(m_tmp_files.back().c_str(), std::ofstream::binary);
	out.write(reinterpret_cast<char const*>(m_buffer.data()), m_buffer.size() * sizeof(T));
	m_buffer.clear();
}

template <typename T>
std::string sorted_external_vector<T>::get_tmp_output_filename(uint64_t id) const {
	std::stringstream filename;
	filename << m_tmp_dirname << "/";
	filename << "lphash.tmp.run";
	if (m_prefix != "") filename << "_" << m_prefix;
	filename << "_" << id << ".bin";
	return filename.str();
}

template <typename T>
sorted_external_vector<T>::const_iterator::const_iterator(sorted_external_vector<T> const* vec) 
    : v(vec), m_mm_files(v->m_tmp_files.size()), heap_idx_comparator([this](uint32_t i, uint32_t j) {return (*m_iterators[i] > *m_iterators[j]);})
{
    m_iterators.reserve(v->m_tmp_files.size());
    m_idx_heap.reserve(v->m_tmp_files.size());

    /* create the input iterators and make the heap */
    for (uint64_t i = 0; i != v->m_tmp_files.size(); ++i) {
        m_mm_files[i].open(v->m_tmp_files.at(i), mm::advice::sequential);
        m_iterators.emplace_back(m_mm_files[i].data(), m_mm_files[i].data() + m_mm_files[i].size());
        m_idx_heap.push_back(i);
    }
    std::make_heap(m_idx_heap.begin(), m_idx_heap.end(), heap_idx_comparator);
}

template <typename T>
sorted_external_vector<T>::const_iterator::const_iterator() : v(nullptr), heap_idx_comparator([]([[maybe_unused]] uint32_t i, [[maybe_unused]] uint32_t j) { return false; }) {}

template <typename T>
T const& sorted_external_vector<T>::const_iterator::operator*() const {
	return *m_iterators[m_idx_heap.front()];
}

template <typename T>
void sorted_external_vector<T>::const_iterator::operator++() {
	advance_heap_head();
}

template <typename T>
bool sorted_external_vector<T>::const_iterator::operator==(const_iterator const& other) const {
	return m_idx_heap.size() == other.m_idx_heap.size(); // TODO make it a little bit stronger
}

template <typename T>
bool sorted_external_vector<T>::const_iterator::operator!=(const_iterator const& other) const {
	return !(operator==(other));
}

template <typename T>
void sorted_external_vector<T>::const_iterator::advance_heap_head() {
	uint32_t idx = m_idx_heap.front();
	++m_iterators[idx];
	if (m_iterators[idx].has_next()) { // percolate down the head
		uint64_t pos = 0;
		uint64_t size = m_idx_heap.size();
		while (2 * pos + 1 < size) {
			uint64_t i = 2 * pos + 1;
			if (i + 1 < size and heap_idx_comparator(m_idx_heap[i], m_idx_heap[i + 1])) ++i;
			if (heap_idx_comparator(m_idx_heap[i], m_idx_heap[pos])) break;
			std::swap(m_idx_heap[pos], m_idx_heap[i]);
			pos = i;
		}
	} else {
		std::pop_heap(m_idx_heap.begin(), m_idx_heap.end(), heap_idx_comparator);
		m_idx_heap.pop_back();
	}
}

} // namespace lphash

#endif