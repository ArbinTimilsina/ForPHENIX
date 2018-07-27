// -*- mode: c++; -*-

#ifndef XJETREC_HISTOGRAM_H_
#define XJETREC_HISTOGRAM_H_

#if defined(__CINT__) || defined(G__DICTIONARY)
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H
#endif // defined(__CINT__) || defined(G__DICTIONARY)

#include <vector>
#include <string>
#include <jetstat/sparse.h>

namespace jet {

	/**
	 *
	 * Histogram with dimension > 2^31 - 1 will lead to undefined
	 * behavior, due to limitations of OpenMP.
	 */
	template<typename x_t, typename c_t>
	class histogram_t {
	public:
		typedef x_t variable_t;
		typedef c_t count_t;
		typedef sparse_hash_map_vector_t<c_t> internal_vector_storage_t;
		typedef sparse_csv_vector_t<c_t> external_vector_storage_t;
		class fixed_bin_t {
		private:
			inline int floor_log2(uint32_t n) const
			{
				int digit = 0;

				if(n >= (1U << 16)) { n >>= 16; digit += 16; }
				if(n >= (1U <<  8)) { n >>=  8; digit +=  8; }
				if(n >= (1U <<  4)) { n >>=  4; digit +=  4; }
				if(n >= (1U <<  2)) { n >>=  2; digit +=  2; }
				if(n >= (1U <<  1)) {           digit +=  1; }

				return n == 0 ? -1 : digit;
			}
		protected:
			std::string _name;
			uint64_t _nbin;
			uint64_t _nbin_with_overflow;
			std::pair<x_t, x_t> _range;
			x_t _width;
			x_t _scale;
		public:
			inline fixed_bin_t(void)
			{
			}
			inline fixed_bin_t(const std::string &name,
							   const uint64_t nbin,
							   const x_t x0, const x_t x1)
				: _name(name), _nbin(nbin),
				  _nbin_with_overflow(nbin + 2), _range(x0, x1)
			{
				_width = _range.second - _range.first;
				_scale = _nbin / _width;
			}
			inline std::string name(void) const
			{
				return _name;
			}
			inline uint64_t nbin(void) const
			{
				return _nbin;
			}
			inline std::pair<x_t, x_t> range(void) const
			{
				return _range;
			}
			inline int nbin_with_overflow(void) const
			{
				return _nbin_with_overflow;
			}
			inline int nbin_with_overflow_padded(void) const
			{
				return floor_log2(_nbin_with_overflow + 1);
			}
			inline x_t scale(void) const
			{
				return _scale;
			}
			inline uint64_t map(const x_t x) const
			{
				if(x < _range.first) {
					return 0;
				}
				else if(x < _range.second) {
					return static_cast<uint64_t>(
						(x - _range.first) * _scale) + 1;
				}
				else {
					return _nbin + 1;
				}
			}
		};
	protected:
		std::string _name;
		std::vector<fixed_bin_t> _fixed_bin;
		internal_vector_storage_t _content;
		internal_vector_storage_t _count;
		uint64_t _ntotal_bin;
		// Member function pointers that provide a primitive form of
		// self-modifying code
		uint64_t (histogram_t<x_t, c_t>::* _map)(
			const std::vector<x_t> &x) const;
		void (histogram_t<x_t, c_t>::* _set_bin_overflow)(
			const c_t content, const c_t count);
		inline void set_bin_overflow_1(const c_t content = 1,
									   const c_t count = 1)
		{
			_content[0] = content;
			_count[0] = count;
			_content[_ntotal_bin - 1] = content;
			_count[_ntotal_bin - 1] = count;
		}
		inline uint64_t map_fixed_0(const std::vector<x_t> &x) const
		{
			return _ntotal_bin - 1;
		}
		inline uint64_t map_fixed_1(const std::vector<x_t> &x) const
		{
			return _fixed_bin[0].map(x[0]);
		}
		inline uint64_t map_fixed_2(const std::vector<x_t> &x) const
		{
			return _fixed_bin[0].map(x[0]) +
				_fixed_bin[0].nbin_with_overflow() *
				_fixed_bin[1].map(x[1]);
		}
		inline uint64_t map_fixed_3(const std::vector<x_t> &x) const
		{
			return _fixed_bin[0].map(x[0]) +
				_fixed_bin[0].nbin_with_overflow() *
				(_fixed_bin[1].map(x[1]) +
				 _fixed_bin[1].nbin_with_overflow() *
				 _fixed_bin[2].map(x[2]));
		}
		inline uint64_t map_fixed_n(const std::vector<x_t> &x) const
		{
			int last = _fixed_bin.size() - 1;
			uint64_t index = _fixed_bin[last].map(x[last]);

			for(int i = last - 1; i >= 0; i--) {
				index *= _fixed_bin[i].nbin_with_overflow();
				index += _fixed_bin[i].map(x[i]);
			}

			return index;
		}
	public:
		inline histogram_t(void)
		{
		}
		inline histogram_t(const std::string &name,
						   const std::vector<fixed_bin_t> &bin)
			: _name(name), _fixed_bin(bin)
		{
			_ntotal_bin = 1;
			for(size_t i = 0; i < _fixed_bin.size(); i++) {
				_ntotal_bin *= _fixed_bin[i].nbin_with_overflow();
			}
			_content = internal_vector_storage_t(_ntotal_bin);
			_count = internal_vector_storage_t(_ntotal_bin);

			switch(_fixed_bin.size()) {
			case 0:
				_map = &histogram_t::map_fixed_0;
				_set_bin_overflow = NULL;
				break;
			case 1:
				_map = &histogram_t::map_fixed_1;
				_set_bin_overflow = &histogram_t::set_bin_overflow_1;
				break;
			case 2:
				_map = &histogram_t::map_fixed_2;
				_set_bin_overflow = NULL;
				break;
			case 3:
				_map = &histogram_t::map_fixed_3;
				_set_bin_overflow = NULL;
				break;
			default:
				_map = &histogram_t::map_fixed_n;
				_set_bin_overflow = NULL;
				break;
			}
		}
		inline histogram_t(const std::string &name,
						   const fixed_bin_t &bin)
			: _name(name),
			  _fixed_bin(std::vector<fixed_bin_t>(1, bin))
		{
			_ntotal_bin = _fixed_bin[0].nbin_with_overflow();
			_content = internal_vector_storage_t(_ntotal_bin);
			_count = internal_vector_storage_t(_ntotal_bin);

			_map = map_fixed_1;
			_set_bin_overflow = set_bin_overflow_1;
		}
		inline std::string name(void) const
		{
			return _name;
		}
		inline std::vector<fixed_bin_t> fixed_bin(void) const
		{
			return _fixed_bin;
		}
		inline internal_vector_storage_t &content(void)
		{
			return _content;
		}
		inline internal_vector_storage_t &count(void)
		{
			return _count;
		}
		inline external_vector_storage_t content_ext(void) const
		{
			return _content;
		}
		inline external_vector_storage_t count_ext(void) const
		{
			return _count;
		}
		inline size_t size(void) const
		{
			return _content.size();
		}
		inline size_t non_zero_size(void) const
		{
			return _content.non_zero_size();
		}
		inline c_t content(const uint64_t n) const
		{
			return _content[n];
		}
		inline c_t count(const uint64_t n) const
		{
			return _count[n];
		}
		inline double error(const uint64_t n) const
		{
			return sqrt(count(n));
		}
		inline std::pair<c_t, c_t> operator[](const uint64_t n) const
		{
			return std::pair<c_t, c_t>(content(n), count(n));
		}
		inline uint64_t flat_index(const std::vector<uint64_t> &n) const
		{
			int last = _fixed_bin.size() - 1;
			uint64_t index = n[last];

			for(int i = last - 1; i >= 0; i--) {
				index *= _fixed_bin[i].nbin_with_overflow();
				index += n[i];
			}

			return index;
		}
		inline std::pair<c_t, c_t>
		operator[](const std::vector<uint64_t> &n) const
		{
			const uint64_t index = flat_index(n);

			return std::pair<c_t, c_t>(content(index), count(index));
		}
		inline void set_bin(const std::vector<uint64_t> &n,
							const c_t content, const c_t count)
		{
			const uint64_t index = flat_index(n);

			_content[index] = content;
			_count[index] = count;
		}
		inline void set_bin_overflow(const c_t content = 1,
									 const c_t count = 1)
		{
			if(_set_bin_overflow == NULL) {
				std::cerr << __FILE__ << ':' << __LINE__
						  << ": error: not implemented" << std::endl;
				return;
			}

			(this->*_set_bin_overflow)(content, count);
		}
		inline void add_to_bin(const std::vector<uint64_t> &n,
							   const c_t content, const c_t count)
		{
			const uint64_t index = flat_index(n);

			_content[index] += content;
			_count[index] += count;
		}
		inline void fill(const std::vector<x_t> &x,
						 const c_t weight = 1)
		{
			uint64_t index = (this->*_map)(x);

			_content[index] += weight;
			_count[index] += weight * weight;
		}
		inline void fill(const x_t x, const c_t weight = 1)
		{
			fill(std::vector<x_t>(1, x), weight);
		}
		inline void fill(const x_t x0, const x_t x1,
						 const c_t weight = 1)
		{
			std::vector<x_t> x;

			x.push_back(x0);
			x.push_back(x1);
			fill(x, weight);
		}
		inline void fill(const x_t x0, const x_t x1, const x_t x2,
						 const c_t weight = 1)
		{
			std::vector<x_t> x;

			x.push_back(x0);
			x.push_back(x1);
			x.push_back(x2);
			fill(x, weight);
		}
		inline histogram_t<x_t, c_t> &
		operator+=(const histogram_t<x_t, c_t> &term)
		{
#ifdef _OPENMP
#pragma omp parallel sections
			{
#pragma omp section
#endif // _OPENMP
				_content += term._content;
#ifdef _OPENMP
#pragma omp section
#endif // _OPENMP
				_count += term._count;
#ifdef _OPENMP
			}
#endif // _OPENMP

			return *this;
		}
		inline histogram_t<x_t, c_t> &
		add_scaled(const histogram_t<x_t, c_t> &x, const c_t a = 1)
		{
#ifdef _OPENMP
#pragma omp parallel sections
			{
#pragma omp section
#endif // _OPENMP
				_content += histogram_t<x_t, c_t>::
					internal_vector_storage_t(x._content).scale(a);
#ifdef _OPENMP
#pragma omp section
#endif // _OPENMP
				_count += histogram_t<x_t, c_t>::
					internal_vector_storage_t(x._count).scale(a * a);
#ifdef _OPENMP
			}
#endif // _OPENMP

			return *this;
		}
		template<typename o_t>
		inline operator histogram_t<x_t, o_t>(void) const
		{
			std::vector<typename histogram_t<x_t, o_t>::
				fixed_bin_t> fixed_bin;

			for(typename std::vector<fixed_bin_t>::const_iterator
					iterator = _fixed_bin.begin();
				iterator != _fixed_bin.end(); iterator++) {
				fixed_bin.push_back(
					histogram_t<x_t, o_t>::fixed_bin_t(
						iterator->name(), iterator->nbin(),
						iterator->range().first,
						iterator->range().second));
			}

			histogram_t<x_t, o_t> ret(_name, fixed_bin);

#ifdef _OPENMP
#pragma omp parallel sections
			{
#pragma omp section
#endif // _OPENMP
				ret.content() = _content;
#ifdef _OPENMP
#pragma omp section
#endif // _OPENMP
				ret.count() = _count;
#ifdef _OPENMP
			}
#endif // _OPENMP

			return ret;
		}
	};

	template<typename c_t>
	inline void fwrite_histogram_array_rle(
		const sparse_csv_vector_t<c_t> &vector, FILE *fp)
	{
		std::vector<uint16_t> repetition(1, 0);
		std::vector<c_t> value(1, vector[0]);

		for(uint64_t i = 1; i < vector.size(); i++) {
			if(vector[i] != value.back() ||
			   repetition.back() == USHRT_MAX) {
				repetition.push_back(0);
				value.push_back(vector[i]);
			}
			else {
				repetition.back()++;
			}
		}

		uint64_t size = repetition.size();

		fwrite(&size, sizeof(uint64_t), 1, fp);
		fwrite(&(repetition[0]), sizeof(uint16_t), size, fp);
		fwrite(&(value[0]), sizeof(c_t), size, fp);
	}

	template<typename c_t>
	inline void fread_histogram_array_rle(
		sparse_hash_map_vector_t<c_t> &vector, FILE *fp)
	{
		uint64_t size;

		fread(&size, sizeof(uint64_t), 1, fp);

		std::vector<uint16_t> repetition(size, 0);
		std::vector<c_t> value(size, vector[0]);

		fread(&(repetition[0]), sizeof(uint16_t), size, fp);
		fread(&(value[0]), sizeof(c_t), size, fp);

		uint64_t vector_size = 0;

		for(std::vector<uint16_t>::const_iterator iterator =
				repetition.begin();
			iterator != repetition.end(); iterator++) {
			vector_size += *iterator + 1;
		}
		vector = sparse_hash_map_vector_t<c_t>(vector_size);

		uint64_t index = 0;

		for(uint64_t i = 0; i < size; i++) {
			if(value[i] != 0) {
				for(uint32_t j = 0; j <= repetition[i]; j++) {
					vector[index] = value[i];
					index++;
				}
			}
			else {
				index += repetition[i] + 1;
			}
		}
	}

	template<typename c_t>
	inline void fwrite_histogram_array_csv(
		const sparse_csv_vector_t<c_t> &vector, FILE *fp)
	{
		std::vector<uint64_t> index;
		std::vector<c_t> value;

		index.reserve(vector.non_zero_size());
		value.reserve(vector.non_zero_size());
		for(typename sparse_csv_vector_t<c_t>::const_iterator
				iterator = vector.begin();
			iterator != vector.end(); iterator++) {
			index.push_back(iterator->first);
			value.push_back(iterator->second);
		}

		uint64_t size[2];

		size[0] = vector.size();
		size[1] = index.size();

		fwrite(size, sizeof(uint64_t), 2, fp);
		fwrite(&(index[0]), sizeof(uint64_t), size[1], fp);
		fwrite(&(value[0]), sizeof(c_t), size[1], fp);
	}

	template<typename c_t>
	inline bool fread_histogram_array_csv(
		sparse_hash_map_vector_t<c_t> &vector, FILE *fp)
	{
		uint64_t size[2];

		if(fread(size, sizeof(uint64_t), 2, fp) != 2) {
			return false;
		}

		std::vector<uint64_t> index(size[1], 0);
		std::vector<c_t> value(size[1], 0);

		if(fread(&(index[0]), sizeof(uint64_t), size[1], fp) != size[1]) {
			return false;
		}
		if(fread(&(value[0]), sizeof(c_t), size[1], fp) != size[1]) {
			return false;
		}

		vector = sparse_hash_map_vector_t<c_t>(size[0]);

		for(uint64_t i = 0; i < size[1]; i++) {
			vector[index[i]] = value[i];
		}

		return true;
	}

#ifdef HAVE_GMM
	template<typename c_t>
	inline void fwrite_histogram_array_rle(
		const gmm::rsvector<c_t> &vector, FILE *fp)
	{
		std::vector<uint16_t> repetition(1, 0);
		std::vector<c_t> value(1, vector[0]);

		for(uint64_t i = 1; i < vector.size(); i++) {
			if(vector[i] != value.back() ||
			   repetition.back() == USHRT_MAX) {
				repetition.push_back(0);
				value.push_back(vector[i]);
			}
			else {
				repetition.back()++;
			}
		}

		uint64_t size = repetition.size();

		fwrite(&size, sizeof(uint64_t), 1, fp);
		fwrite(&(repetition[0]), sizeof(uint16_t), size, fp);
		fwrite(&(value[0]), sizeof(c_t), size, fp);
	}

	template<typename c_t>
	inline void fread_histogram_array_rle(
		gmm::wsvector<c_t> &vector, FILE *fp)
	{
		uint64_t size;

		fread(&size, sizeof(uint64_t), 1, fp);

		std::vector<uint16_t> repetition(size, 0);
		std::vector<c_t> value(size, vector[0]);

		fread(&(repetition[0]), sizeof(uint16_t), size, fp);
		fread(&(value[0]), sizeof(c_t), size, fp);

		uint64_t vector_size = 0;

		for(std::vector<uint16_t>::const_iterator iterator =
				repetition.begin();
			iterator != repetition.end(); iterator++) {
			vector_size += *iterator + 1;
		}
		vector = gmm::wsvector<c_t>(vector_size);

		uint64_t index = 0;

		for(uint64_t i = 0; i < size; i++) {
			if(value[i] != 0) {
				for(uint32_t j = 0; j <= repetition[i]; j++) {
					vector[index] = value[i];
					index++;
				}
			}
			else {
				index += repetition[i] + 1;
			}
		}
	}

	template<typename c_t>
	inline void fwrite_histogram_array_csv(
		const gmm::rsvector<c_t> &vector, FILE *fp)
	{
		std::vector<uint64_t> index;
		std::vector<c_t> value;

		index.reserve(gmm::nnz(vector));
		value.reserve(gmm::nnz(vector));
		for(gmm::rsvector_const_iterator<c_t> iterator =
				vect_const_begin(vector);
			iterator != vect_const_end(vector); iterator++) {
			index.push_back(iterator.index());
			value.push_back(*iterator);
		}

		uint64_t size[2];

		size[0] = vector.size();
		size[1] = index.size();

		fwrite(size, sizeof(uint64_t), 2, fp);
		fwrite(&(index[0]), sizeof(uint64_t), size[1], fp);
		fwrite(&(value[0]), sizeof(c_t), size[1], fp);
	}

	template<typename c_t>
	inline bool fread_histogram_array_csv(
		gmm::wsvector<c_t> &vector, FILE *fp)
	{
		uint64_t size[2];

		if(fread(size, sizeof(uint64_t), 2, fp) != 2) {
			return false;
		}

		std::vector<uint64_t> index(size[1], 0);
		std::vector<c_t> value(size[1], 0);

		if(fread(&(index[0]), sizeof(uint64_t), size[1], fp) != size[1]) {
			return false;
		}
		if(fread(&(value[0]), sizeof(c_t), size[1], fp) != size[1]) {
			return false;
		}

		vector = gmm::wsvector<c_t>(size[0]);

#if 1
		for(uint64_t i = 0; i < size[1]; i++) {
			vector[index[i]] = value[i];
		}
#else
		const uint64_t stride = 20 + 2;

		for(uint64_t i = 0; i < size[1]; i++) {
			uint64_t div = index[i] / stride;
			uint64_t mod = index[i] % stride;

			if(mod >= 1 && mod < 5) {
				const uint64_t reduced_index = div * stride + 1;
				vector[reduced_index] += value[i];
			}
		}
#endif

		return true;
	}
#endif // HAVE_GMM

	template<typename x_t, typename c_t>
	inline void fwrite_histogram_deprecated(
		const histogram_t<x_t, c_t> &histogram, FILE *fp)
	{
		std::vector<typename histogram_t<x_t, c_t>::fixed_bin_t>
			fixed_bin = histogram.fixed_bin();

		uint64_t nfixed_bin = fixed_bin.size();

		fwrite(&nfixed_bin, sizeof(uint64_t), 1, fp);

		for(uint64_t i = 0; i < nfixed_bin; i++) {
			uint64_t nbin = fixed_bin[i].nbin();

			fwrite(&nbin, sizeof(uint64_t), 1, fp);

			x_t range[2];

			range[0] = fixed_bin[i].range().first;
			range[1] = fixed_bin[i].range().second;

			fwrite(range, sizeof(x_t), 2, fp);
		}

		fwrite_histogram_array_rle(histogram.content_ext(), fp);
		fwrite_histogram_array_rle(histogram.count_ext(), fp);
	}

	template<typename x_t, typename c_t>
	inline bool fread_histogram_deprecated(
		histogram_t<x_t, c_t> &histogram, FILE *fp)
	{
		uint64_t nfixed_bin;

		fread(&nfixed_bin, sizeof(uint64_t), 1, fp);

		std::vector<typename histogram_t<x_t, c_t>::fixed_bin_t>
			fixed_bin;

		for(uint64_t i = 0; i < nfixed_bin; i++) {
			uint64_t nbin;
			x_t range[2];

			fread(&nbin, sizeof(uint64_t), 1, fp);
			fread(&range, sizeof(x_t), 2, fp);

			fixed_bin.push_back(
				typename histogram_t<x_t, c_t>::fixed_bin_t(
					"", nbin, range[0], range[1]));
		}

		histogram = histogram_t<x_t, c_t>("", fixed_bin);

		fread_histogram_array_rle(histogram.content(), fp);
		fread_histogram_array_rle(histogram.count(), fp);

		return true;
	}

	template<typename x_t, typename c_t>
	inline void fwrite_histogram(
		const histogram_t<x_t, c_t> &histogram, FILE *fp)
	{
		uint64_t name_size = 0;

		fwrite(&name_size, sizeof(uint64_t), 1, fp);

		std::vector<typename histogram_t<x_t, c_t>::fixed_bin_t>
			fixed_bin = histogram.fixed_bin();

		uint64_t nfixed_bin = fixed_bin.size();

		fwrite(&nfixed_bin, sizeof(uint64_t), 1, fp);

		for(uint64_t i = 0; i < nfixed_bin; i++) {
			fwrite(&name_size, sizeof(uint64_t), 1, fp);

			uint64_t nbin = fixed_bin[i].nbin();

			fwrite(&nbin, sizeof(uint64_t), 1, fp);

			x_t range[2];

			range[0] = fixed_bin[i].range().first;
			range[1] = fixed_bin[i].range().second;

			fwrite(range, sizeof(x_t), 2, fp);
		}

		fwrite_histogram_array_csv(histogram.content_ext(), fp);
		fwrite_histogram_array_csv(histogram.count_ext(), fp);
	}

	template<typename x_t, typename c_t>
	inline bool fread_histogram(
		histogram_t<x_t, c_t> &histogram, FILE *fp)
	{
		uint64_t name_size = 0;

		fread(&name_size, sizeof(uint64_t), 1, fp);

		uint64_t nfixed_bin;

		fread(&nfixed_bin, sizeof(uint64_t), 1, fp);

		std::vector<typename histogram_t<x_t, c_t>::fixed_bin_t>
			fixed_bin;

		for(uint64_t i = 0; i < nfixed_bin; i++) {
			uint64_t nbin;
			x_t range[2];

			fread(&name_size, sizeof(uint64_t), 1, fp);
			fread(&nbin, sizeof(uint64_t), 1, fp);
			fread(&range, sizeof(x_t), 2, fp);

			fixed_bin.push_back(
				typename histogram_t<x_t, c_t>::fixed_bin_t(
					"", nbin, range[0], range[1]));
		}

		histogram = histogram_t<x_t, c_t>("", fixed_bin);

		if(!fread_histogram_array_csv(histogram.content(), fp)) {
			return false;
		}
		if(!fread_histogram_array_csv(histogram.count(), fp)) {
			return false;
		}
		return true;
	}

}

#endif // XJETREC_HISTOGRAM_H_
