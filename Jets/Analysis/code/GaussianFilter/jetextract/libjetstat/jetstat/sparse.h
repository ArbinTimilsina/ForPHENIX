// -*- mode: c++; -*-

#ifndef XJETREC_SPARSE_H_
#define XJETREC_SPARSE_H_

#if defined(__CINT__) || defined(G__DICTIONARY)
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H
#endif // defined(__CINT__) || defined(G__DICTIONARY)

#include <vector>
#include <utility>
//#ifdef HAVE_SPARSEHASH
#include <google/sparse_hash_map>
#include <jetbase/hash.h>
//#endif // HAVE_SPARSEHASH
#include <map>
#ifdef HAVE_GMM
#include <gmm/gmm.h>
#endif // HAVE_GMM
#include <jetbase/dbc.h>

namespace jet {

	// Template sparse histogram class that can handle 1.8 x 10^19
	// elements

	template<typename c_t>
	class sparse_hash_map_vector_t :
		public google::sparse_hash_map<
		uint64_t, c_t, jet::hash<uint64_t> > {
	private:
		inline bool hash_test(void) const
		{
			// This tests the hash function using trivial numbers plus
			// prime numbers close to 2^32, 2^64, and fixed point
			// representations of the inverse golden ratio.

			if(jet::hash<uint64_t>()(0ULL) ==
			   static_cast<size_t>(0x7208f7fa198a2d81ULL) &&
			   jet::hash<uint64_t>()(1ULL) ==
			   static_cast<size_t>(0x8fbb8d815c9e092eULL) &&
			   jet::hash<uint64_t>()(0x9e3779cdULL) ==
			   static_cast<size_t>(0x9d0dc18c5d8bfee5ULL) &&
			   jet::hash<uint64_t>()(0xfffffffbULL) ==
			   static_cast<size_t>(0xc8c599f2981bac4eULL) &&
			   jet::hash<uint64_t>()(0x10000000fULL) ==
			   static_cast<size_t>(0x26fc04b83682bc4eULL) &&
			   jet::hash<uint64_t>()(0x5555555555555555ULL) ==
			   static_cast<size_t>(0xc2127871791aac90ULL) &&
			   jet::hash<uint64_t>()(0x9e3779b97f4a7c55ULL) ==
			   static_cast<size_t>(0xac4f33ac39d660b5ULL) &&
			   jet::hash<uint64_t>()(0xaaaaaaaaaaaaaaaaULL) ==
			   static_cast<size_t>(0x1be750a32849a5feULL) &&
			   jet::hash<uint64_t>()(0xffffffffffffffc5ULL) ==
			   static_cast<size_t>(0x6e5854b9c1f3ad0dULL) &&
			   jet::hash<uint64_t>()(0xffffffffffffffffULL) ==
			   static_cast<size_t>(0x3a8e59c4e2c576ddULL)) {
				return true;
			}
			else {
				std::cerr << __FILE__ << ':' << __LINE__
						  << ": error: unrecognized hash function"
						  << std::endl;
				return false;
			}
		}
	protected:
		uint64_t _dense_size;
	public:
		typedef google::sparse_hash_map<
			uint64_t, c_t, jet::hash<uint64_t> > base_type;
		typedef typename base_type::size_type size_type;
		typedef typename base_type::iterator iterator;
		typedef typename base_type::const_iterator const_iterator;
		inline sparse_hash_map_vector_t(const uint64_t size = 0)
			: _dense_size(size)
		{
			if(!hash_test()) {
				exit(EXIT_FAILURE);
			}
		}
		inline sparse_hash_map_vector_t<c_t>
		operator+=(const sparse_hash_map_vector_t<c_t> &v)
		{
			for(typename base_type::const_iterator iterator =
					v.begin();
				iterator != v.end(); iterator++) {
				(*this)[iterator->first] =
					(*this)[iterator->first] + iterator->second;
			}

			return *this;
		}
		inline uint64_t non_zero_size(void) const
		{
			return base_type::size();
		}
		inline uint64_t size(void) const
		{
			return _dense_size;
		}
		inline sparse_hash_map_vector_t<c_t> scale(const c_t a) const
		{
			sparse_hash_map_vector_t<c_t> ret = *this;

			for(typename base_type::iterator iterator = ret.begin();
				iterator != ret.end(); iterator++) {
				iterator->second *= a;
			}

			return ret;
		}
		template<typename o_t>
		inline operator sparse_hash_map_vector_t<o_t>(void) const
		{
			sparse_hash_map_vector_t<o_t> ret;

			for(typename base_type::const_iterator iterator =
					this->begin();
				iterator != this->end(); iterator++) {
				ret[iterator->first] =
					static_cast<o_t>(iterator->second);
			}

			return ret;
		}
	};

	template<typename c_t>
	class sparse_csv_vector_t :
		public std::map<uint64_t, c_t> {
	protected:
		uint64_t _dense_size;
	public:
		typedef std::map<uint64_t, c_t> base_type;
		typedef typename base_type::size_type size_type;
		typedef typename base_type::iterator iterator;
		typedef typename base_type::const_iterator const_iterator;
		inline sparse_csv_vector_t(const uint64_t size = 0)
			: _dense_size(size)
		{
		}
		inline sparse_csv_vector_t(
			const sparse_hash_map_vector_t<c_t> &v)
			: _dense_size(v.size())
		{
			for(typename sparse_hash_map_vector_t<c_t>::
					const_iterator iterator = v.begin();
				iterator != v.end(); iterator++) {
				(*this)[iterator->first] = iterator->second;
			}
		}
		template<typename vector_t>
		inline sparse_csv_vector_t<c_t>
		operator+=(const vector_t &v)
		{
			for(typename vector_t::const_iterator iterator = v.begin();
				iterator != v.end(); iterator++) {
				(*this)[iterator->first] =
					(*this)[iterator->first] + iterator->second;
			}

			return *this;
		}
		inline uint64_t non_zero_size(void) const
		{
			return base_type::size();
		}
		inline uint64_t size(void) const
		{
			return _dense_size;
		}
		inline operator std::vector<c_t>(void)
		{
			std::vector<c_t> retval(_dense_size, 0);

			for(typename base_type::const_iterator iterator = this->begin();
				iterator != this->end(); iterator++) {
				retval[iterator->first] = iterator->second;
			}

			return retval;
		}
	};

	template<typename v_t>
	class sparse_tensor_t : public v_t {
	protected:
		std::vector<uint64_t> _leading_dimension;
		uint64_t (sparse_tensor_t<v_t>::* _flat_index)(
			const std::vector<uint64_t> &x) const;
		inline uint64_t
		flat_index_0(const std::vector<uint64_t> &index) const
		{
			return 0;
		}
		inline uint64_t
		flat_index_1(const std::vector<uint64_t> &index) const
		{
			return index[0];
		}
		inline uint64_t
		flat_index_2(const std::vector<uint64_t> &index) const
		{
			return index[0] + _leading_dimension[0] *
				index[1];
		}
		inline uint64_t
		flat_index_3(const std::vector<uint64_t> &index) const
		{
			return index[0] + _leading_dimension[0] *
				(index[1] + _leading_dimension[1] *
				 index[2]);
		}
		inline uint64_t
		flat_index_4(const std::vector<uint64_t> &index) const
		{
			return index[0] + _leading_dimension[0] *
				(index[1] + _leading_dimension[1] *
				 (index[2] + _leading_dimension[2] *
				  index[3]));
		}
		inline uint64_t
		flat_index_5(const std::vector<uint64_t> &index) const
		{
			return index[0] + _leading_dimension[0] *
				(index[1] + _leading_dimension[1] *
				 (index[2] + _leading_dimension[2] *
				  (index[3] + _leading_dimension[3] *
				   index[4])));
		}
		inline uint64_t
		flat_index_6(const std::vector<uint64_t> &index) const
		{
			return index[0] + _leading_dimension[0] *
				(index[1] + _leading_dimension[1] *
				 (index[2] + _leading_dimension[2] *
				  (index[3] + _leading_dimension[3] *
				   (index[4] + _leading_dimension[4] *
					index[5]))));
		}
		inline uint64_t
		flat_index_7(const std::vector<uint64_t> &index) const
		{
			return index[0] + _leading_dimension[0] *
				(index[1] + _leading_dimension[1] *
				 (index[2] + _leading_dimension[2] *
				  (index[3] + _leading_dimension[3] *
				   (index[4] + _leading_dimension[4] *
					(index[5] + _leading_dimension[5] *
					 index[6])))));
		}
		inline uint64_t
		flat_index_8(const std::vector<uint64_t> &index) const
		{
			return index[0] + _leading_dimension[0] *
				(index[1] + _leading_dimension[1] *
				 (index[2] + _leading_dimension[2] *
				  (index[3] + _leading_dimension[3] *
				   (index[4] + _leading_dimension[4] *
					(index[5] + _leading_dimension[5] *
					 (index[6] + _leading_dimension[6] *
					  index[6]))))));
		}
		inline uint64_t
		flat_index_9(const std::vector<uint64_t> &index) const
		{
			return index[0] + _leading_dimension[0] *
				(index[1] + _leading_dimension[1] *
				 (index[2] + _leading_dimension[2] *
				  (index[3] + _leading_dimension[3] *
				   (index[4] + _leading_dimension[4] *
					(index[5] + _leading_dimension[5] *
					 (index[6] + _leading_dimension[6] *
					  (index[7] + _leading_dimension[7] *
					   index[8])))))));
		}
		inline uint64_t
		flat_index_10(const std::vector<uint64_t> &index) const
		{
			return index[0] + _leading_dimension[0] *
				(index[1] + _leading_dimension[1] *
				 (index[2] + _leading_dimension[2] *
				  (index[3] + _leading_dimension[3] *
				   (index[4] + _leading_dimension[4] *
					(index[5] + _leading_dimension[5] *
					 (index[6] + _leading_dimension[6] *
					  (index[7] + _leading_dimension[7] *
					   (index[8] + _leading_dimension[8] *
						index[9]))))))));
		}
		inline uint64_t
		flat_index_11(const std::vector<uint64_t> &index) const
		{
			return index[0] + _leading_dimension[0] *
				(index[1] + _leading_dimension[1] *
				 (index[2] + _leading_dimension[2] *
				  (index[3] + _leading_dimension[3] *
				   (index[4] + _leading_dimension[4] *
					(index[5] + _leading_dimension[5] *
					 (index[6] + _leading_dimension[6] *
					  (index[7] + _leading_dimension[7] *
					   (index[8] + _leading_dimension[8] *
						(index[9] + _leading_dimension[9] *
						 index[10])))))))));
		}
		inline uint64_t
		flat_index_12(const std::vector<uint64_t> &index) const
		{
			return index[0] + _leading_dimension[0] *
				(index[1] + _leading_dimension[1] *
				 (index[2] + _leading_dimension[2] *
				  (index[3] + _leading_dimension[3] *
				   (index[4] + _leading_dimension[4] *
					(index[5] + _leading_dimension[5] *
					 (index[6] + _leading_dimension[6] *
					  (index[7] + _leading_dimension[7] *
					   (index[8] + _leading_dimension[8] *
						(index[9] + _leading_dimension[9] *
						 (index[10] + _leading_dimension[10] *
						  index[11]))))))))));
		}
		inline uint64_t
		flat_index_13(const std::vector<uint64_t> &index) const
		{
			return index[0] + _leading_dimension[0] *
				(index[1] + _leading_dimension[1] *
				 (index[2] + _leading_dimension[2] *
				  (index[3] + _leading_dimension[3] *
				   (index[4] + _leading_dimension[4] *
					(index[5] + _leading_dimension[5] *
					 (index[6] + _leading_dimension[6] *
					  (index[7] + _leading_dimension[7] *
					   (index[8] + _leading_dimension[8] *
						(index[9] + _leading_dimension[9] *
						 (index[10] + _leading_dimension[10] *
						  (index[11] + _leading_dimension[11] *
						   index[12])))))))))));
		}
		inline uint64_t
		flat_index_14(const std::vector<uint64_t> &index) const
		{
			return index[0] + _leading_dimension[0] *
				(index[1] + _leading_dimension[1] *
				 (index[2] + _leading_dimension[2] *
				  (index[3] + _leading_dimension[3] *
				   (index[4] + _leading_dimension[4] *
					(index[5] + _leading_dimension[5] *
					 (index[6] + _leading_dimension[6] *
					  (index[7] + _leading_dimension[7] *
					   (index[8] + _leading_dimension[8] *
						(index[9] + _leading_dimension[9] *
						 (index[10] + _leading_dimension[10] *
						  (index[11] + _leading_dimension[11] *
						   (index[12] + _leading_dimension[12] *
							index[13]))))))))))));
		}
		inline uint64_t
		flat_index_15(const std::vector<uint64_t> &index) const
		{
			return index[0] + _leading_dimension[0] *
				(index[1] + _leading_dimension[1] *
				 (index[2] + _leading_dimension[2] *
				  (index[3] + _leading_dimension[3] *
				   (index[4] + _leading_dimension[4] *
					(index[5] + _leading_dimension[5] *
					 (index[6] + _leading_dimension[6] *
					  (index[7] + _leading_dimension[7] *
					   (index[8] + _leading_dimension[8] *
						(index[9] + _leading_dimension[9] *
						 (index[10] + _leading_dimension[10] *
						  (index[11] + _leading_dimension[11] *
						   (index[12] + _leading_dimension[12] *
							(index[13] + _leading_dimension[13] *
							 index[14])))))))))))));
		}
		inline uint64_t
		flat_index_n(const std::vector<uint64_t> &index) const
		{
			uint64_t flat_index = index.back();

			for(size_t i = _leading_dimension.size() - 2;
				i >= 0; i--) {
				flat_index *= _leading_dimension[i];
				flat_index += index[i];
			}

			return flat_index;
		}
		inline void set_execution_path(void)
		{
			switch(_leading_dimension.size()) {
			case 0:
				_flat_index = &sparse_tensor_t::flat_index_0;
				break;
			case 1:
				_flat_index = &sparse_tensor_t::flat_index_1;
				break;
			case 2:
				_flat_index = &sparse_tensor_t::flat_index_2;
				break;
			case 3:
				_flat_index = &sparse_tensor_t::flat_index_3;
				break;
			case 4:
				_flat_index = &sparse_tensor_t::flat_index_4;
				break;
			case 5:
				_flat_index = &sparse_tensor_t::flat_index_5;
				break;
			case 6:
				_flat_index = &sparse_tensor_t::flat_index_6;
				break;
			case 7:
				_flat_index = &sparse_tensor_t::flat_index_7;
				break;
			case 8:
				_flat_index = &sparse_tensor_t::flat_index_8;
				break;
			case 9:
				_flat_index = &sparse_tensor_t::flat_index_9;
				break;
			case 10:
				_flat_index = &sparse_tensor_t::flat_index_10;
				break;
			case 11:
				_flat_index = &sparse_tensor_t::flat_index_11;
				break;
			case 12:
				_flat_index = &sparse_tensor_t::flat_index_12;
				break;
			case 13:
				_flat_index = &sparse_tensor_t::flat_index_13;
				break;
			case 14:
				_flat_index = &sparse_tensor_t::flat_index_14;
				break;
			case 15:
				_flat_index = &sparse_tensor_t::flat_index_15;
				break;
			default:
				_flat_index = &sparse_tensor_t::flat_index_n;
				// _set_bin_overflow = NULL;
				break;
			}
		}
	public:
		inline sparse_tensor_t(void)
			: _leading_dimension(0)
		{
			set_execution_path();
		}
		inline uint64_t
		flat_index(const std::vector<uint64_t> &index) const
		{
			return (this->_flat_index)(index);
		}
		inline typename v_t::base_type::data_type
		operator[](const std::vector<uint64_t> &index) const
		{
			return (*this)[flat_index(index)];
		}
	};

}

#endif // XJETREC_SPARSE_H_
