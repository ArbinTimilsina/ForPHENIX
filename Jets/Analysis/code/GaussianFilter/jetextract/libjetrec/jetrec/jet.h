// -*- mode: c++; -*-

#ifndef XJETREC_JET_H_
#define XJETREC_JET_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <cmath>
#include <climits>
#include <jetevent/particle.h>

#ifndef ULLONG_MAX
#define ULLONG_MAX 18446744073709551615ULL
#endif // ULLONG_MAX

/////////////////////////////////////////////////////////////////////

// Jet definition

namespace jet {

	class reconstruction_filtering_t;

	class jet_t {
	private:
		snowmass_vector_t _momentum;
		std::vector<std::vector<track_t>::const_iterator>
		_constituent;
		snowmass_vector_t _seed_momentum;
		uint64_t _pixel_index;
		float _area;
		float _sigma_perp_per_area;
	public:
		jet_t(void)
			: _momentum(NAN, NAN, NAN, NAN),
			  _seed_momentum(NAN, NAN, NAN, NAN),
			  _pixel_index(ULLONG_MAX), _area(NAN),
			  _sigma_perp_per_area(NAN)
		{
		}
		jet_t(const snowmass_vector_t &momentum)
			: _momentum(momentum), _seed_momentum(momentum),
			  _pixel_index(ULLONG_MAX), _area(NAN),
			  _sigma_perp_per_area(NAN)
		{
		}
		jet_t(const uint64_t pixel_index,
			  const snowmass_vector_t &momentum)
			: _momentum(momentum), _seed_momentum(momentum),
			  _pixel_index(pixel_index), _area(NAN),
			  _sigma_perp_per_area(NAN)
		{
		}
		inline operator snowmass_vector_t(void) const
		{
			return _momentum;
		}
		inline snowmass_vector_t momentum(void) const
		{
			return _momentum;
		}
		inline snowmass_vector_t &momentum(void)
		{
			return _momentum;
		}
		inline std::vector<std::vector<track_t>::const_iterator>
		constituent(void) const
		{
			return _constituent;
		}
		inline std::vector<std::vector<track_t>::const_iterator> &
		constituent(void)
		{
			return _constituent;
		}
		inline snowmass_vector_t seed_momentum(void) const
		{
			return _seed_momentum;
		}
		inline snowmass_vector_t &seed_momentum(void)
		{
			return _seed_momentum;
		}
		inline uint64_t pixel_index(void) const
		{
			return _pixel_index;
		}
		inline uint64_t &pixel_index(void)
		{
			return _pixel_index;
		}
		inline float area(void) const
		{
			return _area;
		}
		inline float &area(void)
		{
			return _area;
		}
		inline float sigma_perp_per_area(void) const
		{
			return _sigma_perp_per_area;
		}
		inline float &sigma_perp_per_area(void)
		{
			return _sigma_perp_per_area;
		}
		inline float radial_seed_distance(void) const
		{
			return _momentum.radial_distance(_seed_momentum);
		}
		inline jet_t recombine(const jet_t &p) const
		{
			jet_t ret;

			ret._momentum = _momentum.recombine(p._momentum);
			ret._constituent.resize(
				_constituent.size() + p._constituent.size());
			std::copy(p._constituent.begin(),
					  p._constituent.end(),
					  std::copy(_constituent.begin(),
								_constituent.end(),
								ret._constituent.begin()));

			return ret;
		}
		inline bool operator<(const jet_t p) const
		{
			return _momentum < p._momentum;
		}
		friend class reconstruction_filtering_t;
		friend class reconstruction_ktjet_t;
		friend class reconstruction_fastjet_t;
		friend class reconstruction_siscone_t;
	};

}

#endif // XJETREC_JET_H_
