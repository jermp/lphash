#include "../include/mphf.hpp"

// #include "../include/prettyprint.hpp"

namespace lphash {

mphf::mphf()
    : k(0)
    , m(0)
    , mm_seed(0)
    , nkmers(0)
    , distinct_minimizers(0)
    , n_maximal(0)
    , right_coll_sizes_start(0)
    , none_sizes_start(0)
    , none_pos_start(0) 
    , max_ram(0)
{
    mphf_configuration.minimal_output = true;
    mphf_configuration.seed = constants::seed;
    mphf_configuration.c = constants::c;
    mphf_configuration.alpha = 0.94;
    mphf_configuration.verbose_output = false;
    mphf_configuration.num_threads = 0;
    mphf_configuration.tmp_dir = "";
    mphf_configuration.ram = 8 * essentials::GB;
};

mphf::mphf(uint8_t klen, uint8_t mm_size, uint64_t seed, uint64_t total_number_of_kmers, double c, 
           uint8_t nthreads, uint8_t max_memory, std::string temporary_directory, bool verbose)
    : k(klen)
    , m(mm_size)
    , mm_seed(seed)
    , nkmers(total_number_of_kmers)
    , distinct_minimizers(0)
    , n_maximal(0)
    , right_coll_sizes_start(0)
    , none_sizes_start(0)
    , none_pos_start(0) 
    , max_ram(max_memory)
{
    mphf_configuration.minimal_output = true;
    mphf_configuration.seed = constants::seed;
    mphf_configuration.c = c;
    mphf_configuration.alpha = 0.94;
    mphf_configuration.verbose_output = verbose;
    mphf_configuration.num_threads = nthreads;
    mphf_configuration.ram = static_cast<uint64_t>(max_memory) * essentials::GB;
    if (temporary_directory != "") {
        mphf_configuration.tmp_dir = temporary_directory;
        essentials::create_directory(temporary_directory);
    }
}

mphf::mm_itr_t::mm_itr_t(sorted_external_vector<mm_triplet_t>::const_iterator& mm_itr) : m_iterator(mm_itr) 
{}

uint64_t mphf::mm_itr_t::operator*() const 
{
    return (*m_iterator).itself;
}

void mphf::mm_itr_t::operator++() 
{
    ++m_iterator;
}

mphf::km_itr_t::km_itr_t(sorted_external_vector<kmer_t>::const_iterator& km_itr) : m_iterator(km_itr) 
{}

kmer_t const& mphf::km_itr_t::operator*() const
{
    return (*m_iterator);
}

void mphf::km_itr_t::operator++()
{
    ++m_iterator;
}

void mphf::build_minimizers_mphf(sorted_external_vector<mm_triplet_t>::const_iterator& mm_itr, std::size_t number_of_minimizers) {
    mm_itr_t dummy_itr(mm_itr);
    minimizer_order.build_in_external_memory(std::move(dummy_itr), number_of_minimizers, mphf_configuration);
}

void mphf::build_fallback_mphf(sorted_external_vector<kmer_t>::const_iterator& km_itr, std::size_t number_of_colliding_kmers) {
    km_itr_t dummy_itr(km_itr);
    fallback_kmer_order.build_in_external_memory(dummy_itr, number_of_colliding_kmers, mphf_configuration);
}

void mphf::build_inverted_index(sorted_external_vector<mm_triplet_t>::const_iterator& start, std::size_t number_of_distinct_minimizers) {
    n_maximal = 0;
    std::size_t colliding_minimizers = 0;
    uint64_t universe = 0;
    quartet_wtree_builder wtb(number_of_distinct_minimizers);
    auto cmp64 = []([[maybe_unused]] uint64_t const& a, [[maybe_unused]] uint64_t const& b) {return false;};
    sorted_external_vector<uint64_t> left_positions(max_ram * essentials::GB / 4, cmp64, mphf_configuration.tmp_dir, get_group_id());
    sorted_external_vector<uint64_t> right_or_collision_sizes(max_ram * essentials::GB / 4, cmp64, mphf_configuration.tmp_dir, get_group_id());
    sorted_external_vector<uint64_t> none_sizes(max_ram * essentials::GB / 4, cmp64, mphf_configuration.tmp_dir, get_group_id());
    sorted_external_vector<uint64_t> none_positions(max_ram * essentials::GB / 4, cmp64, mphf_configuration.tmp_dir, get_group_id());
    for (std::size_t i = 0; i < number_of_distinct_minimizers; ++i) {
        mm_triplet_t mm = *start;
        if (mm.size == 0) {
            wtb.push_back(RIGHT_OR_COLLISION);
            right_or_collision_sizes.push_back(0);
            // universe += 0;
            ++colliding_minimizers;
        } else {
            if (mm.p1 == k - m) {
                if (mm.size == k - m + 1) {
                    wtb.push_back(MAXIMAL);
                    ++n_maximal;
                } else {
                    wtb.push_back(RIGHT_OR_COLLISION);
                    right_or_collision_sizes.push_back(mm.size);
                    universe += mm.size;
                }
            } else {
                if (mm.p1 == mm.size - 1) {
                    wtb.push_back(LEFT);
                    left_positions.push_back(mm.p1 + 1);  // +1 because with p1 == 0 we have 1 k-mer in the prefix sum
                    universe += mm.p1 + 1;
                } else {
                    wtb.push_back(NONE);
                    none_positions.push_back(mm.p1); // here we do not have +1 since p1 != 0 by itself (NONE type of super-k-mer)
                    none_sizes.push_back(mm.size);
                    universe += mm.p1 + mm.size;
                }
            }
        }
        ++start;
    }
    assert(none_positions.size() == none_sizes.size());

    wtree.build(wtb);

    right_coll_sizes_start = left_positions.size();
	none_sizes_start = right_coll_sizes_start + right_or_collision_sizes.size();
	none_pos_start = none_sizes_start + none_sizes.size();

    if (mphf_configuration.verbose_output) {
        double maximal = static_cast<double>(n_maximal) / number_of_distinct_minimizers * 100;
        double left = static_cast<double>(left_positions.size()) / number_of_distinct_minimizers * 100;
        double right = static_cast<double>(right_or_collision_sizes.size() - colliding_minimizers) / number_of_distinct_minimizers * 100;
        double none = static_cast<double>(none_positions.size()) / number_of_distinct_minimizers * 100;
        double ambiguous = static_cast<double>(colliding_minimizers) / number_of_distinct_minimizers * 100;
        std::cerr << "Percentage of maximal super-k-mers: " << maximal << "%\n";
        std::cerr << "Percentage of left-maximal super-k-mers: " << left << "%\n";
        std::cerr << "Percentage of right-maximal super-k-mers : " << right << "%\n";
        std::cerr << "Percentage of unclassified super-k-mers: " << none << "%\n";
        std::cerr << "Percentage of ambiguous minimizers: " << ambiguous << "%\n";
    }

    typedef sorted_external_vector<uint64_t>::const_iterator itr_t;
    std::vector<std::pair<itr_t, itr_t>> sp;
    // The following is needed to build the pairs in-place, since const_iterator is non copy-constructible
    sp.emplace_back(std::piecewise_construct, std::make_tuple(left_positions.cbegin()), std::make_tuple(left_positions.cend()));
    sp.emplace_back(std::piecewise_construct, std::make_tuple(right_or_collision_sizes.cbegin()), std::make_tuple(right_or_collision_sizes.cend()));
    sp.emplace_back(std::piecewise_construct, std::make_tuple(none_sizes.cbegin()), std::make_tuple(none_sizes.cend()));
    sp.emplace_back(std::piecewise_construct, std::make_tuple(none_positions.cbegin()), std::make_tuple(none_positions.cend()));
    append_iterator sp_itr(sp);
    cumulative_iterator c_itr(sp_itr);
    assert(number_of_distinct_minimizers == none_pos_start + n_maximal);
    distinct_minimizers = number_of_distinct_minimizers;
	sizes_and_positions.encode(c_itr, none_pos_start + none_positions.size(), universe);
}

uint64_t mphf::get_minimizer_L0() const noexcept { return distinct_minimizers; }

uint64_t mphf::get_kmer_count() const noexcept { return nkmers; }

uint64_t mphf::num_bits() const noexcept {
	auto mm_mphf_size_bits = minimizer_order.num_bits();
	auto triplet_tree_size_bits = wtree.num_bits();
	auto elias_sequence_size_bits = (sizeof(n_maximal) + sizeof(right_coll_sizes_start) + sizeof(none_sizes_start) + sizeof(none_pos_start)) * 8 + sizes_and_positions.num_bits();
	auto kmer_mphf_size_bits = fallback_kmer_order.num_bits();
	auto total_bit_size = mm_mphf_size_bits + triplet_tree_size_bits + elias_sequence_size_bits + kmer_mphf_size_bits +
						  (sizeof(mphf_configuration) + sizeof(k) + sizeof(m) + sizeof(mm_seed) + sizeof(nkmers) + sizeof(distinct_minimizers)) * 8;
	return total_bit_size;
}

uint64_t mphf::get_minimizer_order(uint64_t mm) const
{
    return minimizer_order(mm);
}

mphf::mm_context_t mphf::query(kmer_t kmer, uint64_t minimizer, uint32_t position) const {
	mm_context_t res;
	uint64_t mp_hash = minimizer_order(minimizer);
    // std::cerr << mp_hash << " ";
	auto [mm_type, mm_type_rank] = wtree.rank_of(mp_hash);
    // std::cerr << mm_type << " " << mm_type_rank << "\n";

    // std::cerr << n_maximal << "\n";

	switch (mm_type) {
		case LEFT:
			res.global_rank = sizes_and_positions.access(mm_type_rank) + (k - m + 1) * n_maximal; // number of left-KMERS before our bucket
			res.local_rank = position;
			res.type = LEFT;
			// std::cerr << "[LEFT] rank = " << mm_type_rank << ", ";
			// std::cerr << "global shift = " << res.global_rank << ", local shift = " << res.local_rank;
			break;
		case RIGHT_OR_COLLISION: {
			auto [val1, val2] = sizes_and_positions.pair(right_coll_sizes_start + mm_type_rank);
            assert(val2 >= val1);
			uint64_t sk_size = val2 - val1;
			if (sk_size == 0) {
				res.global_rank = sizes_and_positions.access(none_pos_start) + (k - m + 1) * n_maximal; // prefix sum of all sizes (sizes of collisions are 0)
				res.local_rank = fallback_kmer_order(kmer);
				res.type = NONE + 1;
				// std::cerr << "[COLLISION] rank = " << none_pos_start << ", global shift = " << res.global_rank << ", "; 
                // std::cerr << "local shift = " << res.local_rank;
			} else {
				res.global_rank = val1 + (k - m + 1) * n_maximal; // global shift
				res.local_rank = k - m - position;				  // local shift
				res.type = RIGHT_OR_COLLISION;					  // in this case it is only RIGHT
				// std::cerr << "[RIGHT] rank = " << right_coll_sizes_start + mm_type_rank << ", ";
				// std::cerr << "global shift = " << res.global_rank << ", ";
				// std::cerr << "local shift = " << res.local_rank;
			}
		} break;
		case MAXIMAL: // easy case
			// maximal k-mer hashes are smaller than those of all the other types
			res.global_rank = (k - m + 1) * mm_type_rank;
			res.local_rank = position;
			res.type = MAXIMAL;
			// std::cerr << "[MAXIMAL] rank = " << mm_type_rank << ", ";
			// std::cerr << "global shift = " << res.global_rank << ", local shift = " << res.local_rank;
			break;
		case NONE: {
			// locpres_hash = sizes_and_positions.access(none_sizes_start + mm_type_rank);  // prefix sum of sizes 
            // locpres_hash += sk_size - position;  // position in the first k-mer - actual position = local shift
			res.global_rank = sizes_and_positions.access(none_sizes_start + mm_type_rank) + (k - m + 1) * n_maximal;
			uint64_t sk_size = sizes_and_positions.diff(none_pos_start + mm_type_rank); // p1 actually
			res.local_rank = sk_size - position;
			res.type = NONE;
			// std::cerr << "[NONE] rank = " << none_sizes_start + mm_type_rank << " = " << none_sizes_start << " + " << mm_type_rank << ", "; 
            // std::cerr << "global rank = " << res.global_rank << ", "; 
            // std::cerr << "local shift = " << res.local_rank << ", p1 = " << sk_size << ", p = " << position;
		} break;
		default: throw std::runtime_error("Unrecognized minimizer type");
	}
	res.hval = res.global_rank + res.local_rank;
    // std::cerr << "; hash value = " << res.hval << "\n";
	return res;
}

void mphf::print_statistics() const noexcept {
	auto mm_mphf_size_bits = minimizer_order.num_bits();
	auto triplet_tree_size_bits = wtree.num_bits();
	auto elias_sequence_size_bits = (sizeof(n_maximal) + sizeof(right_coll_sizes_start) + sizeof(none_sizes_start) + sizeof(none_pos_start)) * 8 + sizes_and_positions.num_bits();
	auto kmer_mphf_size_bits = fallback_kmer_order.num_bits();
	auto total_bit_size = mm_mphf_size_bits + triplet_tree_size_bits + elias_sequence_size_bits + kmer_mphf_size_bits +
						  (sizeof(mphf_configuration) + sizeof(k) + sizeof(m) + sizeof(mm_seed) + sizeof(nkmers) + sizeof(distinct_minimizers)) * 8;
	std::cerr << "Total number of k-mers: " << nkmers << "\n";
	std::cerr << "Minimizer MPHF size in bits : " << mm_mphf_size_bits << " (" << static_cast<double>(mm_mphf_size_bits) / total_bit_size * 100 << "%)\n";
	std::cerr << "\t = " << static_cast<double>(mm_mphf_size_bits) / minimizer_order.num_keys() << " bits/minimizer\n\n";
	std::cerr << "Wavelet tree size in bits : " << triplet_tree_size_bits << " (" << static_cast<double>(triplet_tree_size_bits) / total_bit_size * 100 << "%)\n";
	std::cerr << "\t = " << static_cast<double>(triplet_tree_size_bits) / minimizer_order.num_keys() << " bits/minimizer\n\n";
	std::cerr << "Compressed arrays (EF) : " << elias_sequence_size_bits << " (" << static_cast<double>(elias_sequence_size_bits) / total_bit_size * 100 << "%)\n";
	std::cerr << "\t = " << static_cast<double>(elias_sequence_size_bits) / sizes_and_positions.size() << " bits/offset\n\n";
	std::cerr << "Fallback MPHF : " << kmer_mphf_size_bits << " (" << static_cast<double>(kmer_mphf_size_bits) / total_bit_size * 100 << "%)\n";
	std::cerr << "\t = " << static_cast<double>(kmer_mphf_size_bits) / fallback_kmer_order.num_keys() << " bits/kmer\n\n";
	std::cerr << "Total size in bits : " << total_bit_size << "\n";
	std::cerr << "\tequivalent to : " << static_cast<double>(total_bit_size) / nkmers << " bits/k-mer\n";
	std::cerr << "\n";
}

std::ostream& operator<<(std::ostream& out, mphf const& hf) {
	out << "k = " << static_cast<uint32_t>(hf.k) << "\n";
	out << "m = " << static_cast<uint32_t>(hf.m) << "\n";
	out << "minimizer seed = " << hf.mm_seed << "\n";
	out << "number of k-mers = " << hf.nkmers << "\n";
	out << "distinct minimizers = " << hf.distinct_minimizers << "\n";
	out << "maximal super-k-mers = " << hf.n_maximal << "\n";
	out << "starting index of right and collision sizes = " << hf.right_coll_sizes_start << "\n";
	out << "starting index of none sizes = " << hf.none_sizes_start << "\n";
	out << "starting index of none positions = " << hf.none_pos_start << "\n";
	return out;
}

mphf_alt::mphf_alt() : k(0), m(0), mm_seed(0), nkmers(0), distinct_minimizers(0), max_ram(false) {
    mphf_configuration.minimal_output = true;
    mphf_configuration.seed = constants::seed;
    mphf_configuration.c = constants::c;
    mphf_configuration.alpha = 0.94;
    mphf_configuration.verbose_output = false;
    mphf_configuration.num_threads = 0;
    mphf_configuration.tmp_dir = "";
    mphf_configuration.ram = 8 * essentials::GB;
};

mphf_alt::mphf_alt(uint8_t klen, uint8_t mm_size, uint64_t seed, uint64_t total_number_of_kmers, double c, 
                   uint8_t nthreads, uint8_t max_memory, std::string temporary_directory, bool verbose)
    : k(klen), m(mm_size), mm_seed(seed), nkmers(total_number_of_kmers), distinct_minimizers(0), max_ram(max_memory) {
    mphf_configuration.minimal_output = true;
    mphf_configuration.seed = constants::seed;
    mphf_configuration.c = c;
    mphf_configuration.alpha = 0.94;
    mphf_configuration.verbose_output = verbose;
    mphf_configuration.num_threads = nthreads;
    mphf_configuration.ram = static_cast<uint64_t>(max_ram) * essentials::GB;
    if (temporary_directory != "") {
        mphf_configuration.tmp_dir = temporary_directory;
        essentials::create_directory(temporary_directory);
    }
}

std::unordered_set<uint64_t> mphf_alt::build_index(std::vector<mm_triplet_t>& minimizers) {
    {
        std::vector<uint64_t> random_access_buffer;
        {
            std::unordered_set<uint64_t> unique_minimizers;
            for (auto const& triplet : minimizers) unique_minimizers.insert(triplet.itself);
            distinct_minimizers = unique_minimizers.size();
            for (auto const& mm : unique_minimizers) random_access_buffer.push_back(mm);
        }
        if (max_ram == 0) minimizer_order.build_in_internal_memory(random_access_buffer.begin(), distinct_minimizers, mphf_configuration);
        else minimizer_order.build_in_external_memory(random_access_buffer.begin(), distinct_minimizers, mphf_configuration);
    }
	auto mphf_compare = [this](mm_triplet_t const& a, mm_triplet_t const& b) { return minimizer_order(a.itself) < minimizer_order(b.itself); };
	std::sort(minimizers.begin(), minimizers.end(), mphf_compare);

    // std::vector<uint64_t> colliding_minimizers;
    std::unordered_set<uint64_t> colliding_minimizers;
    std::vector<uint32_t> pos;
    std::vector<uint32_t> siz;
    for (std::size_t i = 0; i < minimizers.size();) {
        std::size_t j;
        for (j = i; j < minimizers.size() && minimizers[i].itself == minimizers[j].itself; ++j) {}
        if (j > (i + 1)) {
            // colliding_minimizers.push_back(minimizers[i].itself);
            colliding_minimizers.insert(minimizers[i].itself);
            pos.push_back(0);
            siz.push_back(0);
        } else {
            pos.push_back(minimizers[i].p1);
            siz.push_back(minimizers[i].size);
        }
        i = j;
    }

	positions.encode(pos.begin(), pos.size());
	sizes.encode(siz.begin(), siz.size());
	num_kmers_in_main_index = sizes.access(sizes.size() - 1);
	return colliding_minimizers;
}

void mphf_alt::build_fallback_mphf(std::vector<kmer_t>& colliding_kmers) {
#ifndef NDEBUG
    std::sort(colliding_kmers.begin(), colliding_kmers.end());
    auto it = std::unique(colliding_kmers.begin(), colliding_kmers.end());
    assert(it == colliding_kmers.end());  // if false: there are some duplicates in unbucketable_kmers
#endif
    if (max_ram) fallback_kmer_order.build_in_internal_memory(colliding_kmers.begin(), colliding_kmers.size(), mphf_configuration);
    else fallback_kmer_order.build_in_external_memory(colliding_kmers.begin(), colliding_kmers.size(), mphf_configuration);
}

uint64_t mphf_alt::get_minimizer_L0() const noexcept { return distinct_minimizers; }

uint64_t mphf_alt::get_kmer_count() const noexcept { return nkmers; }

uint64_t mphf_alt::num_bits() const noexcept {
	auto mm_mphf_size_bits = minimizer_order.num_bits();
	auto positions_size_bits = positions.num_bits();
	auto sizes_size_bits = sizes.num_bits();
	auto kmer_mphf_size_bits = fallback_kmer_order.num_bits();
	auto total_bit_size = mm_mphf_size_bits + positions_size_bits + sizes_size_bits + kmer_mphf_size_bits +
						  (sizeof(mphf_configuration) + sizeof(k) + sizeof(m) + sizeof(mm_seed) + sizeof(nkmers) + sizeof(distinct_minimizers)) * 8;
	return total_bit_size;
}

mphf_alt::mm_context_t mphf_alt::query(kmer_t kmer, uint64_t minimizer, uint32_t position) const {
	mm_context_t res;
	uint64_t index = minimizer_order(minimizer);
	auto [val1, val2] = sizes.pair(index);
	uint64_t size = val2 - val1;
	if (size == 0) {
		res.hval = num_kmers_in_main_index + fallback_kmer_order(kmer);
		res.collision = true;
		return res;
	}
	auto p1 = positions.diff(index);
	res.hval = val1 + p1 - position;
	res.collision = false;
	return res;
}

void mphf_alt::print_statistics() const noexcept {
	auto mm_mphf_size_bits = minimizer_order.num_bits();
	auto positions_size_bits = positions.num_bits();
	auto sizes_size_bits = sizes.num_bits();
	auto kmer_mphf_size_bits = fallback_kmer_order.num_bits();
	auto total_bit_size = mm_mphf_size_bits + positions_size_bits + sizes_size_bits + kmer_mphf_size_bits +
						  (sizeof(mphf_configuration) + sizeof(k) + sizeof(m) + sizeof(mm_seed) + sizeof(nkmers) + sizeof(distinct_minimizers)) * 8;
	std::cerr << "Total number of k-mers: " << nkmers << "\n";
	std::cerr << "Minimizer MPHF size in bits : " << mm_mphf_size_bits << " (" << static_cast<double>(mm_mphf_size_bits) / total_bit_size * 100 << "%)\n";
	std::cerr << "\t = " << static_cast<double>(mm_mphf_size_bits) / minimizer_order.num_keys() << " bits/minimizer\n\n";
	std::cerr << "positions EF sequence size in bits : " << positions_size_bits << " (" << static_cast<double>(positions_size_bits) / total_bit_size * 100 << "%)\n";
	std::cerr << "\t = " << static_cast<double>(positions_size_bits) / minimizer_order.num_keys() << " bits/minimizer\n\n";
	std::cerr << "sizes EF sequence size in bits : " << sizes_size_bits << " (" << static_cast<double>(sizes_size_bits) / total_bit_size * 100 << "%)\n";
	std::cerr << "Fallback MPHF : " << kmer_mphf_size_bits << " (" << static_cast<double>(kmer_mphf_size_bits) / total_bit_size * 100 << "%)\n";
	std::cerr << "\t = " << static_cast<double>(kmer_mphf_size_bits) / fallback_kmer_order.num_keys() << " bits/kmer\n\n";
	std::cerr << "Total size in bits : " << total_bit_size << "\n";
	std::cerr << "\tequivalent to : " << static_cast<double>(total_bit_size) / nkmers << " bits/k-mer\n";
	std::cerr << "\n";
}

std::ostream& operator<<(std::ostream& out, mphf_alt const& hf) {
	out << "k = " << static_cast<uint32_t>(hf.k) << "\n";
	out << "m = " << static_cast<uint32_t>(hf.m) << "\n";
	out << "minimizer seed = " << hf.mm_seed << "\n";
	out << "number of k-mers = " << hf.nkmers << "\n";
	out << "distinct minimizers = " << hf.distinct_minimizers << "\n";
	return out;
}

} // namespace lphash