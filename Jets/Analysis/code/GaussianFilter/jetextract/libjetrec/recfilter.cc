#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H

#include <cstdlib>
#include <functional>
#include <jetrec/rec.h>
#include <jetrec/mem.h>

namespace jet {

	void reconstruction_filtering_t::allocate_host(void)
	{
#ifdef HAVE_SSE
		static const int alignment = 16;
		const int padded_size = (_npixel + 3) & ~3;

		I(padded_size >= _npixel);

		builtin_memalign(
			reinterpret_cast<void **>(&_distribution),
			alignment, padded_size * sizeof(float));
		builtin_memalign(
			reinterpret_cast<void **>(&_intermediate),
			alignment, padded_size * sizeof(float));
		if(_background_perp->factorized()) {
			builtin_memalign(
				reinterpret_cast<void **>(
					&_factorized_background_cache),
				alignment,
				background_model_t::nbin_vertex *
				padded_size * sizeof(float));
			_background_cache = NULL;
		}
		else {
			builtin_memalign(
				reinterpret_cast<void **>(&_background_cache),
				alignment,
				background_model_t::nbin_vertex *
				background_model_t::nbin_centrality *
				background_model_t::nbin_reaction_plane *
				padded_size * sizeof(float));
			_factorized_background_cache = NULL;
		}
		builtin_memalign(
			reinterpret_cast<void **>(&_maximum_2),
			alignment, padded_size * sizeof(float));
		builtin_memalign(
			reinterpret_cast<void **>(&_maximum_3),
			alignment, padded_size * sizeof(float));
		builtin_memalign(
			reinterpret_cast<void **>(&_maximum),
			alignment, padded_size * sizeof(bool));

#else // HAVE_SSE
		_distribution = new float[_npixel];
		_intermediate = new float[_npixel];
		if(_background_perp->factorized()) {
			_factorized_background_cache = new float[
				background_model_t::nbin_vertex * _npixel];
		}
		else {
			_background_cache = new float[
				background_model_t::nbin_vertex *
				background_model_t::nbin_centrality *
				background_model_t::nbin_reaction_plane * _npixel];
		}
		_maximum = new bool[_npixel];
#endif // HAVE_SSE
	}

	void reconstruction_filtering_t::free_host(void)
	{
#ifdef HAVE_SSE
		if(_distribution != NULL) {
			builtin_aligned_free(_distribution);
			_distribution = NULL;
		}
		if(_intermediate != NULL) {
			builtin_aligned_free(_intermediate);
			_intermediate = NULL;
		}
		if(_factorized_background_cache != NULL) {
			builtin_aligned_free(_factorized_background_cache);
			_factorized_background_cache = NULL;
		}
		if(_background_cache != NULL) {
			builtin_aligned_free(_background_cache);
			_background_cache = NULL;
		}
		if(_maximum_2 != NULL) {
			builtin_aligned_free(_maximum_2);
			_maximum_2 = NULL;
		}
		if(_maximum_3 != NULL) {
			builtin_aligned_free(_maximum_3);
			_maximum_3 = NULL;
		}
		if(_maximum != NULL) {
			builtin_aligned_free(_maximum);
			_maximum = NULL;
		}
#else // HAVE_SSE
		if(_distribution != NULL) {
			delete [] _distribution;
			_distribution = NULL;
		}
		if(_intermediate != NULL) {
			delete [] _intermediate;
			_intermediate = NULL;
		}
		if(_factorized_background_cache != NULL) {
			delete [] _factorized_background_cache;
			_factorized_background_cache = NULL;
		}
		if(_background_cache != NULL) {
			delete [] _background_cache;
			_background_cache = NULL;
		}
		if(_maximum != NULL) {
			delete [] _maximum;
			_maximum = NULL;
		}
#endif // HAVE_SSE
	}

	// FIXME: Applying the detailed (i.e. nonflat) background should
	// be implemented here.

	void reconstruction_filtering_t::
	fill(const std::vector<track_t>::const_iterator &track_begin,
		 const std::vector<track_t>::const_iterator &track_end,
		 bool column_occupancy[])
		const
	{
		// The dimensional ordering is:
		//
		// row -> azimuth
		// column -> pseudorapidity
		//
		// The ordering should be chosen to place the higher
		// complexity filtering in rows, which then is done using
		// sparse scalar algorithm. Azimutal filtering is 2 pi + O(u *
		// sigma) as expensive as one unit of pseudorapidity, with
		// typically u ~ 3.
		for(int j = 0; j < _npixel; j++) {
			_distribution[j] = 0;
		}
		for(int j = 0; j < _npixel_pseudorapidity; j++) {
			column_occupancy[j] = false;
		}

		I(_npixel_pseudorapidity <= (1 << 23));
		I(_npixel_azimuth <= (1 << 23));

		const float pseudorapidity_offset =
			_pseudorapidity_range.first;
		const float pseudorapidity_scale =
			(float)_npixel_pseudorapidity /
			(_pseudorapidity_range.second -
			 _pseudorapidity_range.first);
		const float azimuth_offset =
			_azimuth_range.first;
		const float azimuth_scale =
			(float)_npixel_azimuth /
			(_azimuth_range.second - _azimuth_range.first);

		for(std::vector<track_t>::const_iterator track_iterator =
				track_begin;
			track_iterator != track_end; track_iterator++) {
			const float pseudorapidity =
				track_iterator->momentum().pseudorapidity();
			const float raw_azimuth =
				track_iterator->momentum().azimuth();
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 1572)
#endif // __INTEL_COMPILER
			const float azimuth = raw_azimuth == (float)M_PI ?
				(float)(-M_PI) : raw_azimuth;
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif // __INTEL_COMPILER
			const float perp = track_iterator->momentum().perp();

			if(pseudorapidity >= _pseudorapidity_range.first &&
			   pseudorapidity < _pseudorapidity_range.second &&
			   azimuth >= _azimuth_range.first &&
			   azimuth < _azimuth_range.second) {
				const int column_index =
					(int)((pseudorapidity - pseudorapidity_offset) *
						  pseudorapidity_scale);
				const int row_index =
					(int)((azimuth - azimuth_offset) * azimuth_scale);

				if(column_index >= 0 &&
				   column_index < _npixel_pseudorapidity &&
				   row_index >= 0 && row_index < _npixel_azimuth) {
					I(column_index >= 0 &&
					  column_index < _npixel_pseudorapidity);
					I(row_index >= 0 && row_index < _npixel_azimuth);

					const int index =
						_npixel_azimuth * column_index + row_index;

					I(index >= 0 && index < _npixel);

					_distribution[index] += perp;
					column_occupancy[column_index] = true;
				}
			}
		}
	}

	void reconstruction_filtering_t::
	fill(const std::vector<track_t>::const_iterator &track_begin,
		 const std::vector<track_t>::const_iterator &track_end,
		 bool column_occupancy[],
		 const snowmass_vector_t &perturbation) const
	{
		// The dimensional ordering is:
		//
		// row -> azimuth
		// column -> pseudorapidity
		//
		// The ordering should be chosen to place the higher
		// complexity filtering in rows, which then is done using
		// sparse scalar algorithm. Azimutal filtering is 2 pi + O(u *
		// sigma) as expensive as one unit of pseudorapidity, with
		// typically u ~ 3.
		for(int j = 0; j < _npixel; j++) {
			_distribution[j] = 0;
		}
		for(int j = 0; j < _npixel_pseudorapidity; j++) {
			column_occupancy[j] = false;
		}

		I(_npixel_pseudorapidity <= (1 << 23));
		I(_npixel_azimuth <= (1 << 23));

		const float pseudorapidity_offset =
			_pseudorapidity_range.first;
		const float pseudorapidity_scale =
			(float)_npixel_pseudorapidity /
			(_pseudorapidity_range.second -
			 _pseudorapidity_range.first);
		const float azimuth_offset =
			_azimuth_range.first;
		const float azimuth_scale =
			(float)_npixel_azimuth /
			(_azimuth_range.second - _azimuth_range.first);

		for(std::vector<track_t>::const_iterator track_iterator =
				track_begin;
			track_iterator != track_end; track_iterator++) {
			const float pseudorapidity =
				track_iterator->momentum().pseudorapidity() +
				perturbation.pseudorapidity();
			const float raw_azimuth =
				track_iterator->momentum().azimuth() +
				perturbation.azimuth();
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 1572)
#endif // __INTEL_COMPILER
			const float azimuth = raw_azimuth == (float)M_PI ?
				(float)(-M_PI) : raw_azimuth;
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif // __INTEL_COMPILER
			const float perp = track_iterator->momentum().perp();

			if(pseudorapidity >= _pseudorapidity_range.first &&
			   pseudorapidity < _pseudorapidity_range.second &&
			   azimuth >= _azimuth_range.first &&
			   azimuth < _azimuth_range.second) {
				const int column_index =
					(int)((pseudorapidity - pseudorapidity_offset) *
						  pseudorapidity_scale);
				const int row_index =
					(int)((azimuth - azimuth_offset) * azimuth_scale);

				if(column_index >= 0 &&
				   column_index < _npixel_pseudorapidity &&
				   row_index >= 0 && row_index < _npixel_azimuth) {
					I(column_index >= 0 &&
					  column_index < _npixel_pseudorapidity);
					I(row_index >= 0 && row_index < _npixel_azimuth);

					const int index =
						_npixel_azimuth * column_index + row_index;

					I(index >= 0 && index < _npixel);

					_distribution[index] += perp;
					column_occupancy[column_index] = true;
				}
			}
		}
	}

	std::vector<jet_t> reconstruction_filtering_t::
	find_jet_maximum(void)
	{
		const int _stride = _npixel_azimuth;

#ifdef HAVE_SSE
		maximum_map_sse();
#else // HAVE_SSE
		maximum_map();
#endif // HAVE_SSE
		stationary_map();

		const float pseudorapidity_offset =
			_pseudorapidity_range.first;
		const float pseudorapidity_scale =
			(_pseudorapidity_range.second -
			 _pseudorapidity_range.first) /
			(float)_npixel_pseudorapidity;
		const float azimuth_offset =
			_azimuth_range.first;
		const float azimuth_scale =
			(_azimuth_range.second -
			 _azimuth_range.first) /
			(float)_npixel_azimuth;
		std::vector<jet_t> jet;

		static const float epsilon = sqrtf(FLT_EPSILON);
		int index = _stride;

		//jet.clear();
		for(int column_index = 1;
			column_index < _npixel_pseudorapidity - 1;
			column_index++) {
			for(int row_index = 0; row_index < _npixel_azimuth;
				row_index++) {
				// I(index == _stride * column_index + row_index);

				if(_maximum[index] &&
				   _distribution[index] > epsilon) {
					const float perp =
						std::max(0.0F, _distribution[index]);
					const float pseudorapidity =
						column_index * pseudorapidity_scale +
						pseudorapidity_offset;
					const float azimuth =
						row_index * azimuth_scale +
						azimuth_offset;
					snowmass_vector_t momentum =
						jet::snowmass_vector_t(
							perp, perp, pseudorapidity,
							azimuth);

					momentum.set_lightlike_time();
					jet.push_back(momentum);
				}
				index++;
			}
		}

		if(jet.size() > 1) {
			std::sort(jet.begin(), jet.end());
		}

		return jet;
	}

	void reconstruction_filtering_t::background_initialize(void)
	{
		const float pseudorapidity_offset =
			_pseudorapidity_range.first;
		const float pseudorapidity_scale =
			(_pseudorapidity_range.second -
			 _pseudorapidity_range.first) /
			(float)_npixel_pseudorapidity;
		const float azimuth_offset =
			_azimuth_range.first;
		const float azimuth_scale =
			(_azimuth_range.second -
			 _azimuth_range.first) /
			(float)_npixel_azimuth;

		if(_background_perp == NULL) {
			return;
		}
		if(_background_perp->factorized()) {
			for(int vertex_index = 0;
				vertex_index < 10;
				vertex_index++) {
			int index = 0;

			for(int column_index = 0;
				column_index < _npixel_pseudorapidity;
				column_index++) {
				const float pseudorapidity =
					column_index * pseudorapidity_scale +
					pseudorapidity_offset;

				for(int row_index = 0; row_index < _npixel_azimuth;
					row_index++) {
					// I(index == _stride * column_index + row_index);

					const float azimuth =
						row_index * azimuth_scale + azimuth_offset;

					std::vector<double> position;

					position.resize(2);
					position[0] = pseudorapidity;
					position[1] =
						jet::angular_range_reduce(azimuth - 0.5 * M_PI);
					_factorized_background_cache[vertex_index * _npixel + index] =
						static_cast<float>(
										   _background_perp->position_dependence()[vertex_index](
								position));
					index++;
				}
			}
			}
		}
		else {
#if 0
			// FIXME: implement vertex dependence
			for(unsigned int centrality_index = 0;
				centrality_index <
					_background_perp->nbin_centrality;
				centrality_index++) {
				std::cerr << __FILE__ << ':' << __LINE__
						  << ": initialize background cache, "
						  << centrality_index
						  << "% centrality" << std::endl;
				for(unsigned int reaction_plane_index = 0;
					reaction_plane_index <
						_background_perp->nbin_reaction_plane;
					reaction_plane_index++) {
					std::vector<double> position;

					position.resize(2);

					for(int column_index = 0;
						column_index < _npixel_pseudorapidity;
						column_index++) {
						const float pseudorapidity =
							column_index * pseudorapidity_scale +
							pseudorapidity_offset;

						position[0] = pseudorapidity;
						for(int row_index = 0;
							row_index < _npixel_azimuth;
							row_index++) {
							// I(index == _stride * column_index +
							// row_index);

							const float azimuth =
								row_index * azimuth_scale +
								azimuth_offset;

							position[1] = jet::angular_range_reduce(
								azimuth - 0.5 * M_PI);
							_background_cache[index] =
								static_cast<float>(
									_background_perp->
									position_dependence(
										vertex_index,
										centrality_index,
										reaction_plane_index)(
											position));
							index++;
						}
					}
				}
			}
#if 0
			FILE *fp = fopen("data/run_7_au_au_background_cache.bin", "w");

			fwrite(_background_cache, sizeof(float), 101 * 32 * _npixel, fp);
			fclose(fp);
#endif
#endif
		}
	}

	void reconstruction_filtering_t::background_compensate(
		const collision_geometry_t &geometry)
	{
		if(!geometry.in_range()) {
			return;
		}
		if(_background_perp == NULL) {
			return;
		}
		if(_background_perp->factorized()) {
			const float scale =
				geometry.centrality() >= 0.0F &&
				geometry.centrality() <= 100.0F ?
				_background_perp->
				vertex_centrality_dependence(geometry) :
				NAN;
			const unsigned int offset =
				factorized_background_model_t::
				vertex_discretize(geometry) * _npixel;

			for(int index = 0; index < _npixel; index++) {
				_distribution[index] -= scale *
					_factorized_background_cache[offset + index];
			}
		}
		else {
#if 0
			// FIXME
			const unsigned int discrete_centrality =
				static_cast<unsigned int>(rint(centrality));
			unsigned int discrete_reaction_plane_low;
			unsigned int discrete_reaction_plane_high;
			double weight_low;
			double weight_high;

			general_background_model_t::discretize_reaction_plane(
				discrete_reaction_plane_low,
				discrete_reaction_plane_high,
				weight_low, weight_high,
				geometry.reaction_plane());

			const unsigned long index_low =
				discrete_centrality *
				background_model_t::nbin_reaction_plane +
				discrete_reaction_plane_low;
			const unsigned long index_high =
				discrete_centrality *
				background_model_t::nbin_reaction_plane +
				discrete_reaction_plane_high;

			if(!(index_low < 101 * 32 &&
				 index_high < 101 * 32)) {
				return;
			}
			for(int index = 0; index < _npixel; index++) {
				_distribution[index] -=
					(weight_low *
					 _background_cache[index_low * _npixel + index] +
					 weight_high *
					 _background_cache[index_high * _npixel + index]);
			}
#endif
		}
	}

	void reconstruction_filtering_t::
	evaluate_step(
		float &perp, float gradient[], float hessian[],
		const jet_t &jet,
		const std::vector<track_t>::const_iterator &track_begin,
		const std::vector<track_t>::const_iterator &track_end,
		const collision_geometry_t &geometry) const
	{
		switch(_filter_type) {
		case FILTER_TYPE_GAUSSIAN:
#ifdef HAVE_SSE
			evaluate_step_gaussian_sse(
				perp, gradient, hessian, jet,
				track_begin, track_end, geometry);
#else // HAVE_SSE
			evaluate_step_gaussian(
				perp, gradient, hessian, jet,
				track_begin, track_end, geometry);
#endif // HAVE_SSE
			break;
		case FILTER_TYPE_EPANECHNIKOV:
			evaluate_step_epanechnikov(
				perp, gradient, hessian, jet,
				track_begin, track_end, geometry);
			break;
		}
	}

	std::vector<jet_t> reconstruction_filtering_t::
	raw_reconstruct(
		const std::vector<track_t>::const_iterator &track_begin,
		const std::vector<track_t>::const_iterator &track_end,
		const collision_geometry_t &geometry)
	{
		bool column_occupancy[_npixel_pseudorapidity];

		fill(track_begin, track_end, column_occupancy);
		filter(column_occupancy);
		background_compensate(geometry);

		const std::vector<jet_t> jet = find_jet_maximum();

		return jet;
	}

#ifdef __INTEL_COMPILER
#pragma warning(disable: 869)
#endif // __INTEL_COMPILER

	std::vector<jet_t> reconstruction_filtering_t::
	ghost_cut(const std::vector<jet_t> &jet) const
	{
#if 0
		const float scale_pseudorapidity =
			(_pseudorapidity_range.second -
			 _pseudorapidity_range.first) /
			(float)_npixel_pseudorapidity;
		const float scale_azimuth =
			(_azimuth_range.second - _azimuth_range.first) /
			(float)_npixel_azimuth;
		const float scale_max =
			std::max(scale_pseudorapidity, scale_azimuth);
#endif
		std::vector<jet_t> ret;

		for(std::vector<jet_t>::const_iterator iterator_outer =
				jet.begin();
			iterator_outer != jet.end(); iterator_outer++) {
			static const float distance_limit =
				powf(FLT_EPSILON, 2.0F / 3.0F);
			bool isolated = true;

			for(std::vector<jet_t>::const_iterator iterator_inner =
					jet.begin();
				iterator_inner != iterator_outer; iterator_inner++) {
				const float distance =
					iterator_outer->momentum().radial_distance(
						iterator_inner->momentum());

				if(distance < distance_limit) {
#if 0
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 1572)
#endif // __INTEL_COMPILER
					if(distance != 0) {
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif // __INTEL_COMPILER
						std::cerr << __FILE__ << ':' << __LINE__
								  << ": distance = "
								  << mathematica_form(distance)
								  << "; radial_seed_distance = {"
								  << mathematica_form(
										iterator_inner->
										radial_seed_distance())
								  << ','
								  << mathematica_form(
										iterator_outer->
										radial_seed_distance())
								  << "}; perp = {"
								  << mathematica_form(
										iterator_inner->
										seed_momentum().perp())
								  << ','
								  << mathematica_form(
										iterator_outer->
										seed_momentum().perp())
								  << "}; pseudorapidity = {"
								  << mathematica_form(
										iterator_inner->
										seed_momentum().
										pseudorapidity())
								  << ','
								  << mathematica_form(
										iterator_outer->
										seed_momentum().
										pseudorapidity())
								  << "};" << std::endl;
					}
#endif
					isolated = false;
					break;
				}
			}
			if(isolated) {
#if 0
				if(iterator_outer->radial_seed_distance() >
				   2.0F * scale_max) {
					std::cerr << __FILE__ << ':' << __LINE__
							  << ": radial_seed_distance = "
							  << mathematica_form(
									iterator_outer->
									radial_seed_distance())
							  << ';' << std::endl;
				}
#endif
				ret.push_back(*iterator_outer);
			}
		}

		return ret;
	}

	std::vector<jet_t> reconstruction_filtering_t::
	reconstruct(
		const std::vector<track_t>::const_iterator &track_begin,
		const std::vector<track_t>::const_iterator &track_end,
		const collision_geometry_t &geometry,
		const bool refine)
	{
		std::vector<jet_t> raw_jet =
			raw_reconstruct(track_begin, track_end, geometry);

		if(refine) {
			const int nraw_jet = raw_jet.size();
			jet_t *praw_jet = &raw_jet[0];

#ifdef _OPENMP
#pragma omp parallel for shared(praw_jet)
#endif // _OPENMP
			for(int i = 0; i < nraw_jet; i++) {
				refine_jet_maximum(
					praw_jet[i], track_begin, track_end, geometry);
			}
			raw_jet = ghost_cut(raw_jet);
		}

#if 0
		const std::vector<jet_t> jet =
			apply_longitudinal_momentum_fraction_limit(
				raw_jet, track_begin, track_end,
				longitudinal_momentum_fraction_limit);
#endif

		return raw_jet;
	}


	std::vector<jet_t> reconstruction_filtering_t::
	embedded_reconstruct(
		const std::vector<track_t>::const_iterator &track_begin,
		const std::vector<track_t>::const_iterator &track_end,
		const float refine_momentum_perp_limit)
	{
		std::vector<jet_t> raw_jet =
			reconstruct(track_begin, track_end,
						collision_geometry_t(), false);

		const int nraw_jet = raw_jet.size();
		jet_t *praw_jet = &raw_jet[0];

#ifdef _OPENMP
#pragma omp parallel for shared(praw_jet)
#endif // _OPENMP
		for(int i = 0; i < nraw_jet; i++) {
			if(praw_jet[i].momentum().perp() >=
			   refine_momentum_perp_limit) {
				refine_jet_maximum(
					praw_jet[i], track_begin, track_end,
					collision_geometry_t());
			}
		}
		raw_jet = ghost_cut(raw_jet);

		return raw_jet;
	}

#if 1
	void reconstruction_filtering_t::
	longitudinal_momentum_fraction_renormalize(
		float &longitudinal_momentum_fraction,
		const snowmass_vector_t &jet_momentum,
		const snowmass_vector_t &track_momentum) const
	{
		const float factor =
			-0.5F / (_standard_deviation * _standard_deviation);
		const float radial_distance_square =
			track_momentum.radial_distance_square(jet_momentum);

		I(factor < 0);
		I(radial_distance_square >= 0);

		const float scale =
			expf(factor * radial_distance_square);

		I(FRANGE(scale, 0, 1));

		longitudinal_momentum_fraction *= scale;
	}

	std::vector<jet_t> reconstruction_filtering_t::
	apply_longitudinal_momentum_fraction_limit(
		const std::vector<jet_t> raw_jet,
		const std::vector<track_t>::const_iterator &track_begin,
		const std::vector<track_t>::const_iterator &track_end,
		const float longitudinal_momentum_fraction_limit)
		const
	{
		std::vector<jet_t> jet;

		if(raw_jet.empty()) {
			return raw_jet;
		}

		if(longitudinal_momentum_fraction_limit < 0.0F) {
			return jet;
		}
		else if(longitudinal_momentum_fraction_limit >= 1.0F) {
			return raw_jet;
		}

		for(std::vector<jet_t>::const_iterator raw_jet_iterator =
				raw_jet.begin();
			raw_jet_iterator != raw_jet.end(); raw_jet_iterator++) {
			const float magnitude_square =
				raw_jet_iterator->momentum().
				cartesian_magnitude_square();

			float max_longitudinal_momentum_fraction = 0;
			for(std::vector<track_t>::const_iterator track_iterator =
					track_begin;
				track_iterator != track_end; track_iterator++) {
				float longitudinal_momentum_fraction =
					track_iterator->momentum().cartesian_dot(
						raw_jet_iterator->momentum()) /
					magnitude_square;

				IG(FEQ(longitudinal_momentum_fraction,
					   track_iterator->momentum().
					   longitudinal_fraction(*raw_jet_iterator)),
				   std::isfinite(longitudinal_momentum_fraction));

				longitudinal_momentum_fraction_renormalize(
					longitudinal_momentum_fraction,
					raw_jet_iterator->momentum(),
					track_iterator->momentum());
				if(longitudinal_momentum_fraction >
				   max_longitudinal_momentum_fraction)
					max_longitudinal_momentum_fraction =
						longitudinal_momentum_fraction;
			}

			if(max_longitudinal_momentum_fraction <=
			   longitudinal_momentum_fraction_limit)
				jet.push_back(*raw_jet_iterator);
		}

		return jet;
	}
#endif

	void reconstruction_filtering_t::
	filter_spatial_domain(void)
	{
		switch(_filter_type) {
		case FILTER_TYPE_GAUSSIAN:
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": filter type = FILTER_TYPE_GAUSSIAN"
					  << std::endl;
			gaussian_filter_spatial_domain();
			break;
		case FILTER_TYPE_EPANECHNIKOV:
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": filter type = FILTER_TYPE_EPANECHNIKOV"
					  << std::endl;
			epanechnikov_filter_spatial_domain();
			break;
		default:
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": error: unknown filter type "
					  << _filter_type << std::endl;
		}
	}

	void reconstruction_filtering_t::
	gaussian_filter_spatial_domain(void)
	{
		const float discrete_standard_deviation_pseudorapidity =
			_standard_deviation *
			(_npixel_pseudorapidity /
			 (_pseudorapidity_range.second -
			  _pseudorapidity_range.first));
		const float discrete_standard_deviation_azimuth =
			_standard_deviation *
			(_npixel_azimuth /
			 (_azimuth_range.second -
			  _azimuth_range.first));
		const float factor_pseudorapidity =
			-0.5F / (discrete_standard_deviation_pseudorapidity *
					 discrete_standard_deviation_pseudorapidity);
		const float factor_azimuth =
			-0.5F / (discrete_standard_deviation_azimuth *
					 discrete_standard_deviation_azimuth);

		for(int j = 0; j < _npixel_pseudorapidity; j++) {
			const float pseudorapidity =
				(float)(j >= (_npixel_pseudorapidity >> 1) ?
						_npixel_pseudorapidity - j : j);

			for(int i = 0; i < _npixel_azimuth; i++) {
				const int index = j * _npixel_azimuth + i;
				const float azimuth =
					(float)(i >= (_npixel_azimuth >> 1) ?
							_npixel_azimuth - i : i);

				_distribution[index] =
					expf(factor_azimuth * (azimuth * azimuth) +
						 factor_pseudorapidity *
						 (pseudorapidity * pseudorapidity));
			}
		}
	}

	void reconstruction_filtering_t::
	epanechnikov_filter_spatial_domain(void)
	{
		const float radius = 1.4F * _standard_deviation;
		const float discrete_standard_deviation_pseudorapidity =
			radius *
			(_npixel_pseudorapidity /
			 (_pseudorapidity_range.second -
			  _pseudorapidity_range.first));
		const float discrete_standard_deviation_azimuth =
			radius *
			(_npixel_azimuth /
			 (_azimuth_range.second -
			  _azimuth_range.first));
		const float factor_pseudorapidity =
			1.0F / (discrete_standard_deviation_pseudorapidity *
					discrete_standard_deviation_pseudorapidity);
		const float factor_azimuth =
			1.0F / (discrete_standard_deviation_azimuth *
					discrete_standard_deviation_azimuth);

		for(int j = 0; j < _npixel_pseudorapidity; j++) {
			const float pseudorapidity =
				(float)(j >= (_npixel_pseudorapidity >> 1) ?
						_npixel_pseudorapidity - j : j);

			for(int i = 0; i < _npixel_azimuth; i++) {
				const int index = j * _npixel_azimuth + i;
				const float azimuth =
					(float)(i >= (_npixel_azimuth >> 1) ?
							_npixel_azimuth - i : i);

				_distribution[index] = std::max(0.0F, 1.0F -
					(factor_azimuth * (azimuth * azimuth) +
					 factor_pseudorapidity *
					 (pseudorapidity * pseudorapidity)));
			}
		}
	}

	float reconstruction_filtering_t::
	filter_weight(const snowmass_vector_t &fragment,
				  const snowmass_vector_t &jet,
				  const float radial_broadening) const
	{
		switch(_filter_type) {
		case FILTER_TYPE_GAUSSIAN:
			{
				if(std::fpclassify(fragment.perp()) == FP_ZERO ||
				   std::fpclassify(jet.perp()) == FP_ZERO)
					return 0.0F;

				const float scale =  -0.5F /
					(std::fpclassify(radial_broadening) == FP_ZERO ?
					 _standard_deviation * _standard_deviation :
					 _standard_deviation * _standard_deviation +
					 radial_broadening * radial_broadening);
				const float dpseudorapidity =
					jet.pseudorapidity() -
					fragment.pseudorapidity();
				const float dazimuth =
					angular_range_reduce(jet.azimuth() -
										 fragment.azimuth());
				const float distance_square =
					dpseudorapidity * dpseudorapidity +
					dazimuth * dazimuth;

				return expf(scale * distance_square);
			}
		}

		return NAN;
	}

	float reconstruction_filtering_t::
	filter_density(
		const snowmass_vector_t &fragment,
		const std::vector<track_t>::const_iterator &track_begin,
		const std::vector<track_t>::const_iterator &track_end) const
	{
		float filter_perp;
		float __attribute__ ((aligned(8))) gradient[2];
		float __attribute__ ((aligned(16))) hessian[4];

		evaluate_step_gaussian(
			filter_perp, gradient, hessian, fragment,
			track_begin, track_end, collision_geometry_t());

		return filter_perp;
	}

	float reconstruction_filtering_t::
	bayesian_weight(
		const snowmass_vector_t &fragment,
		const snowmass_vector_t &jet,
		const std::vector<track_t>::const_iterator &track_begin,
		const std::vector<track_t>::const_iterator &track_end,
		const float radial_broadening) const
	{
		const float filter_perp =
			filter_density(fragment, track_begin, track_end);
		const float weighted_perp = jet.perp() *
			filter_weight(fragment, jet, radial_broadening);
		return std::fpclassify(filter_perp) == FP_ZERO &&
			std::fpclassify(weighted_perp) == FP_ZERO ? 0.0F :
			weighted_perp / filter_perp;
	}

	lorentz_vector_t reconstruction_filtering_t::
	background_lorentz_jet(
		const snowmass_vector_t &momentum,
		const collision_geometry_t &geometry) const
	{
		const float perp = static_cast<float>((*_background_perp)(
			geometry, momentum.pseudorapidity(),
			momentum.azimuth()));
		const float time = static_cast<float>((*_background_time)(
			geometry, momentum.pseudorapidity(),
			momentum.azimuth()));
		const float z = momentum.pseudorapidity() *
			static_cast<float>((*_background_time)(
				geometry, momentum.pseudorapidity(),
				momentum.azimuth()));
		float x;
		float y;

		sincosf(momentum.azimuth(), &y, &x);
		x *= perp;
		y *= perp;

		return lorentz_vector_t(time, x, y, z);
	}

	snowmass_vector_t reconstruction_filtering_t::
	bayesian_jet(
		const snowmass_vector_t &jet,
		const std::vector<track_t>::const_iterator &track_begin,
		const std::vector<track_t>::const_iterator &track_end,
		const collision_geometry_t &geometry,
		const float transverse_momentum_broadening) const
	{
		const float radial_broadening =
			std::fpclassify(transverse_momentum_broadening) ==
			FP_ZERO || std::fpclassify(jet.perp()) == FP_ZERO ?
			0.0F : transverse_momentum_broadening / jet.perp();
		lorentz_vector_t sum(0, 0, 0, 0);

		for(std::vector<track_t>::const_iterator iterator =
				track_begin;
			iterator != track_end; iterator++) {
			if(iterator->is_final_state()) {
				sum += lorentz_vector_t(iterator->momentum()) *
					bayesian_weight(
						iterator->momentum(), jet,
						track_begin, track_end, radial_broadening);
			}
		}

		const lorentz_vector_t ret = sum -
			background_lorentz_jet(jet, geometry);

		return ret;
	}

	lorentz_vector_t reconstruction_filtering_t::
	lorentz_jet(
		const snowmass_vector_t &jet,
		const std::vector<track_t>::const_iterator &track_begin,
		const std::vector<track_t>::const_iterator &track_end,
		const collision_geometry_t &geometry,
		const float transverse_momentum_broadening) const
	{
		lorentz_vector_t sum(0, 0, 0, 0);

		for(std::vector<track_t>::const_iterator iterator =
				track_begin;
			iterator != track_end; iterator++) {
			if(iterator->is_final_state()) {
				sum += lorentz_vector_t(iterator->momentum()) *
					filter_weight(
						iterator->momentum(), jet,
						transverse_momentum_broadening);
			}
		}

		const lorentz_vector_t ret = sum -
			background_lorentz_jet(jet, geometry);

		return ret;
	}

	snowmass_vector_t reconstruction_filtering_t::
	first_order_jet(
		const snowmass_vector_t &jet,
		const std::vector<track_t>::const_iterator &track_begin,
		const std::vector<track_t>::const_iterator &track_end,
		const float transverse_momentum_broadening) const
	{
		float perp = 0;

		for(std::vector<track_t>::const_iterator iterator =
				track_begin;
			iterator != track_end; iterator++) {
			if(iterator->is_final_state()) {
				const snowmass_vector_t fragment =
					iterator->momentum();

				perp += fragment.perp() *
					fragment.radial_distance_square(jet) *
					filter_weight(
						fragment, jet,
						transverse_momentum_broadening);
			}
		}

		snowmass_vector_t ret = jet;

		ret.perp() = perp;
		ret.set_lightlike_time();

		return ret;
	}

	void reconstruction_filtering_t::
	orthogonal_statistics(float &mean, float &rms,
						  const float azimuth) const
	{
		static const float pi_1_3 = (float)(M_PI / 3.0);
		static const float pi_2_3 = (float)(M_PI * 2.0 / 3.0);

		float sum_1 = 0;
		float sum_2 = 0;
		int count = 0;

		const int iazimuth_plus_first = static_cast<int>(
			_npixel_azimuth *
			(angular_range_reduce(azimuth + pi_1_3) -
			 _azimuth_range.first) /
			(_azimuth_range.second - _azimuth_range.first));
		const int iazimuth_plus_second = static_cast<int>(
			_npixel_azimuth *
			(angular_range_reduce(azimuth + pi_2_3) -
			 _azimuth_range.first) /
			(_azimuth_range.second - _azimuth_range.first));

		if(iazimuth_plus_second > iazimuth_plus_first) {
			for(int i = iazimuth_plus_first;
				i < iazimuth_plus_second; i++)
				for(int j = 0; j < _npixel_pseudorapidity; j++) {
					const float u = _distribution[
						j * _npixel_pseudorapidity + i];

					sum_1 += u;
					sum_2 += u * u;
					count++;
				}
		}
		else {
			for(int i = iazimuth_plus_first; i < _npixel_azimuth;
				i++)
				for(int j = 0; j < _npixel_pseudorapidity; j++) {
					const float u = _distribution[
						j * _npixel_pseudorapidity + i];

					sum_1 += u;
					sum_2 += u * u;
					count++;
				}
			for(int i = 0; i < iazimuth_plus_second; i++)
				for(int j = 0; j < _npixel_pseudorapidity; j++) {
					const float u = _distribution[
						j * _npixel_pseudorapidity + i];

					sum_1 += u;
					sum_2 += u * u;
					count++;
				}
		}

		const int iazimuth_minus_first = static_cast<int>(
			_npixel_azimuth *
			(angular_range_reduce(azimuth - pi_2_3) -
			 _azimuth_range.first) /
			(_azimuth_range.second - _azimuth_range.first));
		const int iazimuth_minus_second = static_cast<int>(
			_npixel_azimuth *
			(angular_range_reduce(azimuth - pi_1_3) -
			 _azimuth_range.first) /
			(_azimuth_range.second - _azimuth_range.first));

		if(iazimuth_minus_second > iazimuth_minus_first) {
			for(int i = iazimuth_minus_first;
				i < iazimuth_minus_second; i++)
				for(int j = 0; j < _npixel_pseudorapidity; j++) {
					const float u = _distribution[
						j * _npixel_pseudorapidity + i];

					sum_1 += u;
					sum_2 += u * u;
					count++;
				}
		}
		else {
			for(int i = iazimuth_minus_first;
				i < _npixel_azimuth; i++)
				for(int j = 0; j < _npixel_pseudorapidity; j++) {
					const float u = _distribution[
						j * _npixel_pseudorapidity + i];

					sum_1 += u;
					sum_2 += u * u;
					count++;
				}
			for(int i = 0; i < iazimuth_minus_second; i++)
				for(int j = 0; j < _npixel_pseudorapidity; j++) {
					const float u = _distribution[
						j * _npixel_pseudorapidity + i];

					sum_1 += u;
					sum_2 += u * u;
					count++;
				}
		}

		fprintf(stderr, "count = %d\n", count);
		mean = sum_1 / count;
		rms = count == 0 ? 0 : sqrtf(sum_2 / count);
	}

	double pole_standard_deviation(
		const double q,
		const std::vector<double_complex_t> & pole)
	{
		const int n = pole.size();
		const double q_scale = q / 2;
		double_complex_t cs = 0.0 + 0.0 * _Complex_I;

		for (int i = 0; i < n; i++) {
			const double a = pow(cabs(pole[i]), -1 / q_scale);
			const double t = carg(pole[i]) / q_scale;
			const double_complex_t b = a * cexp(_Complex_I * t);
			const double_complex_t one_b = 1.0 - b;
			const double_complex_t one_b_2 = one_b * one_b;

			cs += 2.0 * b / one_b_2;
		}

		return sqrt(creal(cs));
	}

	double pole_standard_deviation_target(
		const double q,
		const std::vector<double_complex_t> &coefficient)
	{
		std::vector<double_complex_t> pole = coefficient;

		pole.pop_back();

		const double standard_deviation =
			creal(*(coefficient.rbegin()));

		return pole_standard_deviation(q, pole) - standard_deviation;
	}

	const std::string reconstruction_filtering_iir_t::_revision =
		std::string("20081101-IIR");

	std::vector<double_complex_t>
	reconstruction_filtering_iir_t::
	adapt_pole(const double standard_deviation,
			   const std::vector<double_complex_t> &pole)
	{
		// Upper and lower bounds, valid for all standard_deviation >=
		// 0
		const double bound_low = standard_deviation <= 0.2 ?
			0 : 0.8997808 * (standard_deviation - 0.2);
		const double bound_high =
			0.8997810 * (standard_deviation + 0.9);

		std::vector<double_complex_t> coefficient = pole;

		coefficient.push_back(standard_deviation + 0.0 * _Complex_I);

#if 0
		// This can be used to extract the bounds above
		static const int npoint = (1 << 16);
		static const std::pair<double, double> q_range =
			std::pair<double, double>(0.0, 128.0);
		FILE *fp = fopen("output/standard_deviation.txt", "w");

		I(fp != NULL);

		for(int i = 0; i < npoint; i++) {
			const double q = q_range.first +
				(double)i * ((q_range.second - q_range.first) /
							 npoint);
			const double s = pole_standard_deviation(q, pole);

			if(std::isfinite(s)) {
				fprintf(fp, "%.48f %.18f\n", s, q);
			}
		}
		fclose(fp);
#endif

		const double q =
			jet::solve_1d_aps(pole_standard_deviation_target,
							  coefficient, bound_low, bound_high);

		std::cerr << __FILE__ << ':' << __LINE__
				  << ": information: discrete standard deviation = "
				  << standard_deviation << ", q = " << q
				  << std::endl;

		// Adjust poles.
		const int n = pole.size();
		std::vector<double_complex_t> pole_adapted;

		pole_adapted.resize(n);
		for(int i = 0; i < n; i++) {
			const double a = pow(cabs(pole[i]), 2.0 / q);
			const double t = carg(pole[i]) * 2.0 / q;

			pole_adapted[i] = 1.0 / (a * cexp(_Complex_I * t));
		}

		return pole_adapted;
	}

	double reconstruction_filtering_iir_t::
	pole_gain(const std::vector<double_complex_t> &pole)
	{
		const int n = pole.size();
		double_complex_t cg = 1.0 + 0.0 * _Complex_I;

		for(int i = 0; i < n; i++) {
			cg *= 1.0 - pole[i];
		}

		return creal(cg);
	}

	// Calculate the denominator factorization

	double_complex_t reconstruction_filtering_iir_t::
	gr(const double_complex_t pole_j,
	   const std::vector<double_complex_t> &pole,
	   double gain)
	{
		const int n = pole.size();
		double_complex_t conj_pole_j = conj(pole_j);
		double_complex_t inverse_pole_j = 1.0 / pole_j;
		double_complex_t gp = 1.0 + 0.0 * _Complex_I;

		for(int i = 0; i < n; i++) {
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 1572)
#endif // __INTEL_COMPILER
			if(pole[i] != pole_j && pole[i] != conj_pole_j) {
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif // __INTEL_COMPILER
				gp *= 1.0 - pole[i] * inverse_pole_j;
			}
			gp *= 1.0 - pole[i] * pole_j;
		}

		return gain / gp;
	}

	void reconstruction_filtering_iir_t::
	initialize_filter(float coefficient[],
					  const double standard_deviation,
					  const std::vector<double_complex_t> &initial_pole)
	{
		// _g = new Recursive2ndOrderFilter[3][2][2];

		// Adjust the poles for the scale factor q
		const std::vector<double_complex_t> pole =
			adapt_pole(standard_deviation, initial_pole);

		// Filter gain
		const double gain = pole_gain(pole);
		const double gain_2 = gain * gain;

		// Residues
		double_complex_t d0 = pole[0];
		double_complex_t d1 = pole[2];
		double_complex_t e0 = 1.0 / d0;
		double_complex_t e1 = 1.0 / d1;
		double_complex_t g0 = gr(d0, pole, gain_2);
		double_complex_t g1 = gr(d1, pole, gain_2);

		// Coefficients for 2nd-order recursive filters
		double a10 = -2 * creal(d0);
		double a11 = -2 * creal(d1);
		double a20 = cabs(d0) * cabs(d0);
		double a21 = cabs(d1) * cabs(d1);
		double b00;
		double b01;
		double b10;
		double b11;
		double b20;
		double b21;

		// 0th- and 2nd-derivative filters are symmetric
		b10 = cimag(g0) / cimag(e0);
		b11 = cimag(g1) / cimag(e1);
		b00 = creal(g0) - b10 * creal(e0);
		b01 = creal(g1) - b11 * creal(e1);
		b20 = 0;
		b21 = 0;

		// First filter, applied forward
		coefficient[0] = static_cast<float>(b00);
		coefficient[1] = static_cast<float>(b10);
		coefficient[2] = static_cast<float>(b20);
		coefficient[3] = static_cast<float>(a10);
		coefficient[4] = static_cast<float>(a20);

		// Third filter, applied forward
		coefficient[10] = static_cast<float>(b01);
		coefficient[11] = static_cast<float>(b11);
		coefficient[12] = static_cast<float>(b21);
		coefficient[13] = static_cast<float>(a11);
		coefficient[14] = static_cast<float>(a21);

		b20 -= b00 * a20;
		b21 -= b01 * a21;
		b10 -= b00 * a10;
		b11 -= b01 * a11;
		b00 = 0;
		b01 = 0;

		// Second filter, applied backward
		coefficient[5] = static_cast<float>(b00);
		coefficient[6] = static_cast<float>(b10);
		coefficient[7] = static_cast<float>(b20);
		coefficient[8] = static_cast<float>(a10);
		coefficient[9] = static_cast<float>(a20);

		// Fourth filter, applied backward
		coefficient[15] = static_cast<float>(b01);
		coefficient[16] = static_cast<float>(b11);
		coefficient[17] = static_cast<float>(b21);
		coefficient[18] = static_cast<float>(a11);
		coefficient[19] = static_cast<float>(a21);
	}

	void reconstruction_filtering_iir_t::
	filter(const bool column_occupancy[])
	{
		const int _stride = _npixel_azimuth;
		float *pdistribution = _distribution;
		float *pintermediate = _intermediate;

		// Sparse filtering in column (pseudorapidity) is more
		// effective than in row (azimuth).

#ifdef _OPENMP
#pragma omp parallel for shared(pdistribution, pintermediate)
#endif // _OPENMP
		for(int i = 0; i < _npixel_pseudorapidity; i++) {
			if(column_occupancy[i]) {
#ifdef HAVE_SSE
				apply_forward_aos_2(
					pintermediate + _stride * i,
					pdistribution + _stride * i,
					_npixel_azimuth,
					_coefficient_azimuth,
					_coefficient_azimuth + 10);
				accumulate_backward_aos_2(
					pintermediate + _stride * i,
					pdistribution + _stride * i,
					_npixel_azimuth,
					_coefficient_azimuth + 5,
					_coefficient_azimuth + 15);
#else // HAVE_SSE
				apply_forward(
					pintermediate + _stride * i,
					pdistribution + _stride * i,
					_npixel_azimuth,
					_coefficient_azimuth,
					1, true);
				accumulate_backward(
					pintermediate + _stride * i,
					pdistribution + _stride * i,
					_npixel_azimuth,
					_coefficient_azimuth + 5,
					1, true);
				accumulate_forward(
					pintermediate + _stride * i,
					pdistribution + _stride * i,
					_npixel_azimuth,
					_coefficient_azimuth + 10,
					1, true);
				accumulate_backward(
					pintermediate + _stride * i,
					pdistribution + _stride * i,
					_npixel_azimuth,
					_coefficient_azimuth + 15,
					1, true);
#endif // HAVE_SSE
			}
			else {
				for(int j = 0; j < _npixel_azimuth; j++) {
					pintermediate[i * _stride + j] = 0;
				}
			}
		}

#ifdef _OPENMP
#pragma omp parallel for shared(pdistribution, pintermediate)
#endif // _OPENMP
#ifdef HAVE_SSE
		for(int i = 0; i < _npixel_azimuth; i += 4) {
			apply_forward_soa_4(
				pdistribution + i,
				pintermediate + i,
				_npixel_pseudorapidity,
				_coefficient_pseudorapidity,
				_stride);
			accumulate_backward_soa_4(
				pdistribution + i,
				pintermediate + i,
				_npixel_pseudorapidity,
				_coefficient_pseudorapidity + 5,
				_stride);
			accumulate_forward_soa_4(
				pdistribution + i,
				pintermediate + i,
				_npixel_pseudorapidity,
				_coefficient_pseudorapidity + 10,
				_stride);
			accumulate_backward_soa_4(
				pdistribution + i,
				pintermediate + i,
				_npixel_pseudorapidity,
				_coefficient_pseudorapidity + 15,
				_stride);
		}
#else // HAVE_SSE
		for(int i = 0; i < _npixel_azimuth; i++) {
			apply_forward(
				pdistribution + i,
				pintermediate + i,
				_npixel_pseudorapidity,
				_coefficient_pseudorapidity,
				_stride);
			accumulate_backward(
				pdistribution + i,
				pintermediate + i,
				_npixel_pseudorapidity,
				_coefficient_pseudorapidity + 5,
				_stride);
			accumulate_forward(
				pdistribution + i,
				pintermediate + i,
				_npixel_pseudorapidity,
				_coefficient_pseudorapidity + 10,
				_stride);
			accumulate_backward(
				pdistribution + i,
				pintermediate + i,
				_npixel_pseudorapidity,
				_coefficient_pseudorapidity + 15,
				_stride);
		}
#endif // HAVE_SSE

		for(int j = 0; j < _npixel; j++) {
			pdistribution[j] *= _normalization;
		}
	}

	reconstruction_filtering_iir_t::
	reconstruction_filtering_iir_t(
		const float standard_deviation,
		const std::pair<float, float> pseudorapidity_range,
		const int npixel_pseudorapidity,
		const background_model_t &background_perp,
		const background_model_t &background_time,
		const background_model_t &background_z,
		const float *background_cache,
		const int npixel_azimuth,
		const std::pair<float, float> azimuth_range)
		: reconstruction_filtering_t(
			standard_deviation, pseudorapidity_range,
			npixel_pseudorapidity,
			background_perp, background_time, background_z,
			background_cache, npixel_azimuth, azimuth_range)
	{
		_filter_type = FILTER_TYPE_GAUSSIAN;

		static const double_complex_t _anticausal_pole[2] = {
			1.1178503510108764 + 1.2806487439746588 * _Complex_I,
			1.7737391243526701 + 0.4478723498713723 * _Complex_I
		};
		std::vector<double_complex_t> initial_pole;

		initial_pole.push_back(_anticausal_pole[0]);
		initial_pole.push_back(conj(_anticausal_pole[0]));
		initial_pole.push_back(_anticausal_pole[1]);
		initial_pole.push_back(conj(_anticausal_pole[1]));

		const double discrete_standard_deviation_azimuth =
			standard_deviation *
			(_npixel_azimuth /
			 (_azimuth_range.second -
			  _azimuth_range.first));
		const double discrete_standard_deviation_pseudorapidity =
			standard_deviation *
			(_npixel_pseudorapidity /
			 (_pseudorapidity_range.second -
			  _pseudorapidity_range.first));

		_normalization = static_cast<float>(
			2 * M_PI *
			discrete_standard_deviation_pseudorapidity *
			discrete_standard_deviation_azimuth);
		initialize_filter(
			_coefficient_azimuth,
			discrete_standard_deviation_azimuth,
			initial_pole);
		initialize_filter(
			_coefficient_pseudorapidity,
			discrete_standard_deviation_pseudorapidity,
			initial_pole);
	}

	const std::string reconstruction_filtering_fir_t::_revision =
		std::string("20070829-FIR");

#ifdef HAVE_MKL
	void reconstruction_filtering_fir_t::
	dfti_permuted_format_multiply(float to[], const float from[],
								  const int m, const int n) const
	{
		// Note that the MKL uses a column-major ordering
		// (Fortran-like)

		// The size of the upper left corner in the (m, n) matrix is
		// (2 - (n & 1), 2 - (m & 1))

		const int m_2 = (m << 1);
		const int m_n = m * n;

		to[0] *= from[0];
		if((n & 1) == 0) {
			// n % 2 == 0
			to[m] *= from[m];
			for(int index = m_2; index < m_n; index += m_2)
				dfti_packed_multiply(to, from, index, m);
			if((m & 1) == 0) {
				// m % 2 == 0, n % 2 == 0
				to[1] *= from[1];
				to[m + 1] *= from[m + 1];
				for(int index = m_2 + 1; index < m_n; index += m_2)
					dfti_packed_multiply(to, from, index, m);
			}
		}
		else {
			// n % 2 == 1
			for(int index = m; index < m_n; index += m_2)
				dfti_packed_multiply(to, from, index, m);
			if((m & 1) == 0) {
				// m % 2 == 0, n % 2 == 1
				to[1] *= from[1];
				for(int index = m + 1; index < m_n; index += m_2)
					dfti_packed_multiply(to, from, index, m);
			}
		}

		// The contigeous block
		const int row_offset = 2 - (m & 1);

		IG(row_offset == 2, m % 2 == 0);
		IG(row_offset == 1, m % 2 == 1);

		int index = row_offset;

		for(int column = 0; column < n; column++) {
			for(int row = row_offset; row < m; row += 2) {
				I(index >= 0 && index < m_n);
				I(index == column * m + row);

				dfti_packed_multiply(to, from, index, 1);
				index += 2;
			}
			index += row_offset;
		}
	}

#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 869)
#endif // __INTEL_COMPILER
	void reconstruction_filtering_fir_t::
	dfti_filter(const bool column_occupancy[])
	{
		long status;

		status = DftiComputeForward(_dfti_descriptor_handle,
									_distribution, _intermediate);
		if(DftiErrorClass(status, DFTI_NO_ERROR) == 0)
			dfti_error(status);
		dfti_permuted_format_multiply(_intermediate,
									  _filter_fourier_domain,
									  _npixel_azimuth,
									  _npixel_pseudorapidity);
		status = DftiComputeBackward(_dfti_descriptor_handle,
									 _intermediate, _distribution);
		if(DftiErrorClass(status, DFTI_NO_ERROR) == 0)
			dfti_error(status);
	}
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif // __INTEL_COMPILER

	void reconstruction_filtering_fir_t::
	initialize_dfti_filter(void)
	{
		// DFTI_LENGTHS = {m, n}
		long length[2] = {
			_npixel_pseudorapidity, _npixel_azimuth
		};
		long status;

		// Create a DFTI descriptor
		status = DftiCreateDescriptor(&_dfti_descriptor_handle,
									  DFTI_SINGLE, DFTI_REAL,
									  2, length);
		if(DftiErrorClass(status, DFTI_NO_ERROR) == 0)
			dfti_error(status);

		// The DFTI_PERM_FORMAT is likely to be more efficient

		// Set output to be packed complex conjugate-symmetric
		status = DftiSetValue(_dfti_descriptor_handle,
							  DFTI_PACKED_FORMAT, DFTI_PERM_FORMAT);
		if(DftiErrorClass(status, DFTI_NO_ERROR) == 0)
			dfti_error(status);
		// Set placement to be out of place
		status = DftiSetValue(_dfti_descriptor_handle,
							  DFTI_PLACEMENT, DFTI_NOT_INPLACE);
		if(DftiErrorClass(status, DFTI_NO_ERROR) == 0)
			dfti_error(status);

		// DFTI_INPUT_STRIDES = {first_index, stride_in_m,
		// stride_in_n}
		long strides_in[3] = {0, _npixel_azimuth, 1};
		long strides_out[3] = {0, _npixel_azimuth, 1};

		// Set input strides
		status = DftiSetValue(_dfti_descriptor_handle,
							  DFTI_INPUT_STRIDES, strides_in);
		if(DftiErrorClass(status, DFTI_NO_ERROR) == 0)
			dfti_error(status);
		// Set output strides
		status = DftiSetValue(_dfti_descriptor_handle,
							  DFTI_OUTPUT_STRIDES, strides_out);
		if(DftiErrorClass(status, DFTI_NO_ERROR) == 0)
			dfti_error(status);

		const float scale = 1.0F / _npixel;

		status = DftiSetValue(_dfti_descriptor_handle,
							  DFTI_BACKWARD_SCALE, scale);
		if(DftiErrorClass(status, DFTI_NO_ERROR) == 0)
			dfti_error(status);
		// Commit the descriptor
		DftiCommitDescriptor(_dfti_descriptor_handle);
		if(DftiErrorClass(status, DFTI_NO_ERROR) == 0)
			dfti_error(status);

		_filter_fourier_domain = new float[_npixel];

		filter_spatial_domain();
		status = DftiComputeForward(_dfti_descriptor_handle,
									_distribution,
									_filter_fourier_domain);
		if(DftiErrorClass(status, DFTI_NO_ERROR) == 0)
			dfti_error(status);
	}
#endif // HAVE_MKL

#ifdef HAVE_CUDA
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 869)
#endif // __INTEL_COMPILER
	void reconstruction_filtering_fir_t::
	cufft_filter(const bool column_occupancy[])
	{
		cuda_memcpy(_device_distribution, _distribution,
					_npixel * sizeof(cufftReal),
					cudaMemcpyHostToDevice);
		cufft_exec_r2c(_cufft_plan[0], _device_distribution,
					   _device_intermediate);
		cufft_exec_c2r(_cufft_plan[0], _device_intermediate,
					   _device_distribution);
		cuda_memcpy(_distribution, _device_distribution,
					_npixel * sizeof(cufftReal),
					cudaMemcpyDeviceToHost);
	}
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif // __INTEL_COMPILER

	void reconstruction_filtering_fir_t::
	initialize_cufft_filter(void)
	{
		cuda_malloc(reinterpret_cast<void **>
					(&_device_distribution),
					_npixel * sizeof(cufftReal));
		cuda_malloc(reinterpret_cast<void **>
					(&_device_intermediate),
					((_npixel_azimuth >> 1) + 1) *
					((_npixel_pseudorapidity >> 1) + 1) *
					sizeof(cufftComplex));
		cuda_malloc(reinterpret_cast<void **>
					(&_device_filter_fourier_domain),
					((_npixel_azimuth >> 1) + 1) *
					((_npixel_pseudorapidity >> 1) + 1) *
					sizeof(cufftComplex));
		cufft_plan_2d(&_cufft_plan[0], _npixel_azimuth,
					  _npixel_pseudorapidity, CUFFT_R2C);
		filter_spatial_domain();
		cuda_memcpy(_device_distribution, _distribution,
					_npixel * sizeof(cufftReal),
					cudaMemcpyHostToDevice);
		cufft_exec_r2c(_cufft_plan[0], _device_distribution,
					   _device_filter_fourier_domain);
		cufft_plan_2d(&_cufft_plan[1], _npixel_azimuth,
					  _npixel_pseudorapidity, CUFFT_C2R);
	}
#endif // HAVE_CUDA

}
