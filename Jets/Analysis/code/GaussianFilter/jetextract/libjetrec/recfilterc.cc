#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H

#include <cstdlib>
#include <jetrec/rec.h>

namespace jet {

	void reconstruction_filtering_t::
	maximum_map(void) const
	{
		const int _stride = _npixel_azimuth;

		// The naive local maximum finder seems to be faster on the
		// x87 FPU.

		// FIXME: I need to unroll this and see how fast it really can
		// be.

		int index = _stride;

		for(int column_index = 1;
			column_index < _npixel_pseudorapidity - 1;
			column_index++) {
			for(int row_index = 0; row_index < _npixel_azimuth;
				row_index++) {
				I(index == _stride * column_index + row_index);

				const float d =
					_distribution[index];

				const int index_left = row_index == 0 ?
					index + _npixel_azimuth - 1 : index - 1;
				const int index_right =
					row_index == _npixel_azimuth - 1 ?
					index - _npixel_azimuth + 1 : index + 1;
				const int index_top = index - _npixel_azimuth;
				const int index_bottom = index + _npixel_azimuth;
				const int index_top_left = row_index == 0 ?
					index - 1 : index - _npixel_azimuth - 1;
				const int index_top_right =
					row_index == _npixel_azimuth - 1 ?
					index - _npixel_azimuth * 2 + 1 :
					index - _npixel_azimuth + 1;
				const int index_bottom_left = row_index == 0 ?
					index + _npixel_azimuth * 2 - 1 :
					index + _npixel_azimuth - 1;
				const int index_bottom_right =
					row_index == _npixel_azimuth - 1 ?
					index + 1 : index + _npixel_azimuth + 1;

				_maximum[index] =
					d > _distribution[index_top_left] &&
					d > _distribution[index_top] &&
					d > _distribution[index_top_right] &&
					d > _distribution[index_left] &&
					d > _distribution[index_right] &&
					d > _distribution[index_bottom_left] &&
					d > _distribution[index_bottom] &&
					d > _distribution[index_bottom_right];

				// I(_maximum[index] == _maximum_ref[index]);
				index++;
			}
		}
	}

#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 1572)
#endif // __INTEL_COMPILER
	void reconstruction_filtering_t::
	stationary_map(void) const
	{
		const int _stride = _npixel_azimuth;

#if 0
		for(int column_index = 0;
			column_index < _npixel_pseudorapidity;
			column_index++) {
			for(int row_index = 0; row_index < _npixel_azimuth;
				row_index++) {
				const int index =
					_stride * column_index + row_index;
				_maximum[index] = 0;
			}
		}
#endif

		for(int column_index = 1;
			column_index < _npixel_pseudorapidity;
			column_index++) {
			for(int row_index = 0; row_index < _npixel_azimuth;
				row_index++) {
				const int index_0 =
					_stride * column_index + row_index;
				const int index_1 =
					_stride * (column_index - 1) + row_index;

				if(_distribution[index_0] == _distribution[index_1] &&
				   _distribution[index_0] > 0) {
#if 0
					std::cerr << __FILE__ << ':' << __LINE__
							  << ": warning: stationary point at {"
							  << column_index << ',' << row_index
							  << "}, {" << column_index - 1 << ','
							  << row_index << '}' << std::endl;
#endif
					_maximum[index_0] = true;
					_maximum[index_1] = true;
				}
			}
		}
		for(int column_index = 0;
			column_index < _npixel_pseudorapidity;
			column_index++) {
			if(_distribution[_stride * column_index] ==
			   _distribution[_stride * column_index +
							 _npixel_azimuth - 1] &&
			   _distribution[_stride * column_index] > 0) {
#if 0
				std::cerr << __FILE__ << ':' << __LINE__
						  << ": warning: stationary point at {"
						  << column_index << ',' << 0 << "}, {"
						  << column_index << ','
						  << _npixel_azimuth - 1 << '}' << std::endl;
#endif
				_maximum[_stride * column_index] = true;
				_maximum[_stride * column_index +
						 _npixel_azimuth - 1] = true;
			}
			for(int row_index = 1; row_index < _npixel_azimuth;
				row_index++) {
				const int index_0 =
					_stride * column_index + row_index;
				const int index_1 =
					_stride * column_index + row_index - 1;

				if(_distribution[index_0] == _distribution[index_1] &&
				   _distribution[index_0] > 0) {
#if 0
					std::cerr << __FILE__ << ':' << __LINE__
							  << ": warning: stationary point at {"
							  << column_index << ',' << row_index
							  << "}, {" << column_index << ','
							  << row_index - 1 << '}' << std::endl;
#endif
					_maximum[index_0] = true;
					_maximum[index_1] = true;
				}
			}
		}
	}
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif // __INTEL_COMPILER

	void reconstruction_filtering_t::
	evaluate_step_gaussian(
		float &perp, float gradient[], float hessian[],
		const jet_t &jet,
		const std::vector<track_t>::const_iterator &track_begin,
		const std::vector<track_t>::const_iterator &track_end,
		const collision_geometry_t &geometry) const
	{
		const float scale = -0.5F /
			(_standard_deviation * _standard_deviation);
		const int size = track_end - track_begin;
		// OpenMP accumulation variables (which cannot be
		// pointers/arrays)
		float mp_perp = 0;
		float mp_gradient_0 = 0;
		float mp_gradient_1 = 0;
		float mp_hessian_0 = 0;
		float mp_hessian_1 = 0;
		float mp_hessian_2 = 0;

#ifdef _OPENMP
#pragma omp parallel for								\
	reduction(+: mp_perp, mp_gradient_0, mp_gradient_1,	\
			  mp_hessian_0, mp_hessian_1, mp_hessian_2)
#endif // _OPENMP
		// OpenMP requires (for Fortran compatibility) that the loop
		// variable is of type int (and certainly no C++ iterators)
		for(int i = 0; i < size; i++) {
			const float dpseudorapidity =
				jet.momentum().pseudorapidity() -
				track_begin[i].momentum().pseudorapidity();
			const float dazimuth =
				angular_range_reduce(
					jet.momentum().azimuth() -
					track_begin[i].momentum().azimuth());
			const float distance_square =
				dpseudorapidity * dpseudorapidity +
				dazimuth * dazimuth;
			const float exponential_factor =
				track_begin[i].momentum().perp() *
				expf(scale * distance_square);

			mp_perp += exponential_factor;

			const float two_scale = 2.0F * scale;
			const float two_scale_dazimuth = two_scale * dazimuth;
			const float two_scale_dpseudorapidity =
				two_scale * dpseudorapidity;

			mp_gradient_0 +=
				two_scale_dpseudorapidity * exponential_factor;
			mp_gradient_1 +=
				two_scale_dazimuth * exponential_factor;
			mp_hessian_0 += two_scale *
				(1.0F + two_scale_dpseudorapidity * dpseudorapidity) *
				exponential_factor;
			mp_hessian_2 +=
				two_scale_dazimuth * two_scale_dpseudorapidity *
				exponential_factor;
			mp_hessian_1 += two_scale *
				(1.0F + two_scale_dazimuth * dazimuth) *
				exponential_factor;
		}

		perp = mp_perp;
		gradient[0] = mp_gradient_0;
		gradient[1] = mp_gradient_1;
		hessian[0] = mp_hessian_0;
		hessian[1] = mp_hessian_1;
		hessian[2] = mp_hessian_2;

		double perp_background;
		double gradient_background[2];
		double hessian_background[3];

		_background_perp->evaluate_order_2(
			perp_background, gradient_background, hessian_background,
			geometry, jet.momentum().pseudorapidity(),
			jet.momentum().azimuth());

		perp -= (float)perp_background;
		gradient[0] -= (float)gradient_background[0];
		gradient[1] -= (float)gradient_background[1];
		hessian[0] -= (float)hessian_background[0];
		hessian[1] -= (float)hessian_background[1];
		hessian[2] -= (float)hessian_background[2];
	}

	void reconstruction_filtering_t::
	evaluate_step_epanechnikov(
		float &perp, float gradient[], float hessian[],
		const jet_t &jet,
		const std::vector<track_t>::const_iterator &track_begin,
		const std::vector<track_t>::const_iterator &track_end,
		const collision_geometry_t &geometry) const
	{
		const float radius = 1.4F * _standard_deviation;
		const float radius_square = radius * radius;
		const float scale = -1.0F / radius_square;
		const int size = track_end - track_begin;
		// OpenMP accumulation variables (which cannot be
		// pointers/arrays)
		float mp_perp = 0;
		float mp_gradient_0 = 0;
		float mp_gradient_1 = 0;
		float mp_hessian_0 = 0;
		float mp_hessian_1 = 0;

#ifdef _OPENMP
#pragma omp parallel for								\
	reduction(+: mp_perp, mp_gradient_0, mp_gradient_1,	\
			  mp_hessian_0, mp_hessian_1)
#endif // _OPENMP
		// OpenMP requires (for Fortran compatibility) that the loop
		// variable is of type int (and certainly no C++ iterators)
		for(int i = 0; i < size; i++) {
			const float dpseudorapidity =
				jet.momentum().pseudorapidity() -
				track_begin[i].momentum().pseudorapidity();
			const float dazimuth =
				angular_range_reduce(
					jet.momentum().azimuth() -
					track_begin[i].momentum().azimuth());
			const float distance_square =
				dpseudorapidity * dpseudorapidity +
				dazimuth * dazimuth;

			if(distance_square < radius_square) {
				const float momentum_perp =
					track_begin[i].momentum().perp();
				const float two_perp_scale =
					2.0F * momentum_perp * scale;

				mp_perp += momentum_perp *
					(1.0F + scale * distance_square);
				mp_gradient_0 += two_perp_scale * dpseudorapidity;
				mp_gradient_1 += two_perp_scale * dazimuth;
				mp_hessian_0 += two_perp_scale;
				mp_hessian_1 += two_perp_scale;
			}
		}

		perp = mp_perp;
		gradient[0] = mp_gradient_0;
		gradient[1] = mp_gradient_1;
		hessian[0] = mp_hessian_0;
		hessian[1] = mp_hessian_1;
		hessian[2] = 0;

		double perp_background;
		double gradient_background[2];
		double hessian_background[3];

		_background_perp->evaluate_order_2(
			perp_background, gradient_background, hessian_background,
			geometry, jet.momentum().pseudorapidity(),
			jet.momentum().azimuth());

		perp -= (float)perp_background;
		gradient[0] -= (float)gradient_background[0];
		gradient[1] -= (float)gradient_background[1];
		hessian[0] -= (float)hessian_background[0];
		hessian[1] -= (float)hessian_background[1];
		hessian[2] -= (float)hessian_background[2];
	}

	void reconstruction_filtering_iir_t::
	apply_forward(float output[], const float input[], const int n,
				  const float *coefficient, const int stride,
				  const bool cyclic)
	{
		I(stride >= 1);
		I(input != NULL);
		I(output != NULL);

		float yim2 = 0.0F;
		float yim1 = 0.0F;
		float xim2 = 0.0F;
		float xim1 = 0.0F;

		if(stride == 1) {
			if(cyclic) {
#if 0
				for(int i = 0; i < n; i++) {
					const float xi = input[i];
					const float yi = coefficient[0] * xi +
						coefficient[1] * xim1 +
						coefficient[2] * xim2 -
						coefficient[3] * yim1 -
						coefficient[4] * yim2;

					yim2 = yim1;
					yim1 = yi;
					xim2 = xim1;
					xim1 = xi;
				}
#else
				if(std::fpclassify(coefficient[0]) == FP_ZERO)
					for(int i = 0; i < n; i += 2) {
						float yi;

						yi =
							coefficient[1] * xim1 +
							coefficient[2] * xim2 -
							coefficient[3] * yim1 -
							coefficient[4] * yim2;

						yim2 = yi;
						xim2 = input[i];

						yi =
							coefficient[1] * xim2 +
							coefficient[2] * xim1 -
							coefficient[3] * yim2 -
							coefficient[4] * yim1;

						yim1 = yi;
						xim1 = input[i + 1];
					}
				else
					for(int i = 0; i < n; i += 2) {
						float xi;
						float yi;

						xi = input[i];
						yi =
							coefficient[0] * xi +
							coefficient[1] * xim1 -
							coefficient[3] * yim1 -
							coefficient[4] * yim2;

						yim2 = yi;
						xim2 = xi;

						xi = input[i + 1];
						yi =
							coefficient[0] * xi +
							coefficient[1] * xim2 -
							coefficient[3] * yim2 -
							coefficient[4] * yim1;

						yim1 = yi;
						xim1 = xi;
					}
#endif
			}
#if 0
			for(int i = 0; i < n; i++) {
				const float xi = input[i];
				const float yi = coefficient[0] * xi +
					coefficient[1] * xim1 +
					coefficient[2] * xim2 -
					coefficient[3] * yim1 -
					coefficient[4] * yim2;

				output[i] = yi;
				yim2 = yim1;
				yim1 = yi;
				xim2 = xim1;
				xim1 = xi;
			}
#else
			if(std::fpclassify(coefficient[0]) == FP_ZERO)
				for(int i = 0; i < n; i += 2) {
					float yi;

					yi =
						coefficient[1] * xim1 +
						coefficient[2] * xim2 -
						coefficient[3] * yim1 -
						coefficient[4] * yim2;

					output[i] = yi;
					yim2 = yi;
					xim2 = input[i];

					yi =
						coefficient[1] * xim2 +
						coefficient[2] * xim1 -
						coefficient[3] * yim2 -
						coefficient[4] * yim1;

					output[i + 1] = yi;
					yim1 = yi;
					xim1 = input[i + 1];
				}
			else
				for(int i = 0; i < n; i += 2) {
					float xi;
					float yi;

					xi = input[i];
					yi =
						coefficient[0] * xi +
						coefficient[1] * xim1 -
						coefficient[3] * yim1 -
						coefficient[4] * yim2;

					output[i] = yi;
					yim2 = yi;
					xim2 = xi;

					xi = input[i + 1];
					yi =
						coefficient[0] * xi +
						coefficient[1] * xim2 -
						coefficient[3] * yim2 -
						coefficient[4] * yim1;

					output[i + 1] = yi;
					yim1 = yi;
					xim1 = xi;
				}
#endif
		}
		else {
			I(cyclic == false);

			int i_stride = 0;

			for(int i = 0; i < n; i++) {
				const float xi = input[i_stride];
				const float yi = coefficient[0] * xi +
					coefficient[1] * xim1 +
					coefficient[2] * xim2 -
					coefficient[3] * yim1 -
					coefficient[4] * yim2;

				output[i_stride] = yi;
				i_stride += stride;
				yim2 = yim1;
				yim1 = yi;
				xim2 = xim1;
				xim1 = xi;
			}
		}
	}

	void reconstruction_filtering_iir_t::
	accumulate_forward(float output[], const float input[],
					   const int n, const float *coefficient,
					   const int stride, const bool cyclic)
	{
		I(stride >= 1);
		I(input != NULL);
		I(output != NULL);

		float yim2 = 0.0F;
		float yim1 = 0.0F;
		float xim2 = 0.0F;
		float xim1 = 0.0F;

		if(stride == 1) {
			if(cyclic) {
#if 0
				for(int i = 0; i < n; i++) {
					const float xi = input[i];
					const float yi = coefficient[0] * xi +
						coefficient[1] * xim1 +
						coefficient[2] * xim2 -
						coefficient[3] * yim1 -
						coefficient[4] * yim2;

					yim2 = yim1;
					yim1 = yi;
					xim2 = xim1;
					xim1 = xi;
				}
#else
				if(std::fpclassify(coefficient[0]) == FP_ZERO)
					for(int i = 0; i < n; i += 2) {
						float yi;

						yi =
							coefficient[1] * xim1 +
							coefficient[2] * xim2 -
							coefficient[3] * yim1 -
							coefficient[4] * yim2;

						yim2 = yi;
						xim2 = input[i];

						yi =
							coefficient[1] * xim2 +
							coefficient[2] * xim1 -
							coefficient[3] * yim2 -
							coefficient[4] * yim1;

						yim1 = yi;
						xim1 = input[i + 1];
					}
				else
					for(int i = 0; i < n; i += 2) {
						float xi;
						float yi;

						xi = input[i];
						yi =
							coefficient[0] * xi +
							coefficient[1] * xim1 -
							coefficient[3] * yim1 -
							coefficient[4] * yim2;

						yim2 = yi;
						xim2 = xi;

						xi = input[i + 1];
						yi =
							coefficient[0] * xi +
							coefficient[1] * xim2 -
							coefficient[3] * yim2 -
							coefficient[4] * yim1;

						yim1 = yi;
						xim1 = xi;
					}
#endif
			}
			for(int i = 0; i < n; i++) {
				const float xi = input[i];
				const float yi = coefficient[0] * xi +
					coefficient[1] * xim1 +
					coefficient[2] * xim2 -
					coefficient[3] * yim1 -
					coefficient[4] * yim2;

				output[i] += yi;
				yim2 = yim1;
				yim1 = yi;
				xim2 = xim1;
				xim1 = xi;
			}
		}
		else {
			I(cyclic == false);

			int i_stride = 0;

			for(int i = 0; i < n; i++) {
				const float xi = input[i_stride];
				const float yi = coefficient[0] * xi +
					coefficient[1] * xim1 +
					coefficient[2] * xim2 -
					coefficient[3] * yim1 -
					coefficient[4] * yim2;

				output[i_stride] += yi;
				i_stride += stride;
				yim2 = yim1;
				yim1 = yi;
				xim2 = xim1;
				xim1 = xi;
			}
		}
	}

	void reconstruction_filtering_iir_t::
	accumulate_backward(float output[], const float input[],
						const int n, const float *coefficient,
						const int stride, const bool cyclic)
	{
		I(stride >= 1);
		I(input != NULL);
		I(output != NULL);

		float xip1 = 0.0F;
		float xip2 = 0.0F;
		float yip1 = 0.0F;
		float yip2 = 0.0F;

		if(stride == 1) {
			if(cyclic) {
#if 0
				for(int i = n - 1; i >= 0; i--) {
					const float xi = input[i];
					const float yi = coefficient[0] * xi +
						coefficient[1] * xip1 +
						coefficient[2] * xip2 -
						coefficient[3] * yip1 -
						coefficient[4] * yip2;

					yip2 = yip1;
					yip1 = yi;
					xip2 = xip1;
					xip1 = xi;
				}
#else
				if(std::fpclassify(coefficient[0]) == FP_ZERO)
					for(int i = n - 1; i >= 0; i -= 2) {
						float yi;

						yi =
							coefficient[1] * xip1 +
							coefficient[2] * xip2 -
							coefficient[3] * yip1 -
							coefficient[4] * yip2;

						yip2 = yi;
						xip2 = input[i];

						yi =
							coefficient[1] * xip2 +
							coefficient[2] * xip1 -
							coefficient[3] * yip2 -
							coefficient[4] * yip1;

						yip1 = yi;
						xip1 = input[i - 1];
					}
				else
					for(int i = n - 1; i >= 0; i -= 2) {
						float xi;
						float yi;

						xi = input[i];
						yi =
							coefficient[0] * xi +
							coefficient[1] * xip1 -
							coefficient[3] * yip1 -
							coefficient[4] * yip2;

						yip2 = yi;
						xip2 = xi;

						xi = input[i - 1];
						yi =
							coefficient[0] * xi +
							coefficient[1] * xip2 -
							coefficient[3] * yip2 -
							coefficient[4] * yip1;

						yip1 = yi;
						xip1 = xi;
					}
#endif
			}
			for(int i = n - 1; i >= 0; i--) {
				const float xi = input[i];
				const float yi = coefficient[0] * xi +
					coefficient[1] * xip1 +
					coefficient[2] * xip2 -
					coefficient[3] * yip1 -
					coefficient[4] * yip2;

				output[i] += yi;
				yip2 = yip1;
				yip1 = yi;
				xip2 = xip1;
				xip1 = xi;
			}
		}
		else {
			I(cyclic == false);

			int i_stride = (n - 1) * stride;

			for(int i = n - 1; i >= 0; i--) {
				const float xi = input[i_stride];
				const float yi = coefficient[0] * xi +
					coefficient[1] * xip1 +
					coefficient[2] * xip2 -
					coefficient[3] * yip1 -
					coefficient[4] * yip2;

				output[i_stride] += yi;
				i_stride -= stride;
				yip2 = yip1;
				yip1 = yi;
				xip2 = xip1;
				xip1 = xi;
			}
		}
	}

}
