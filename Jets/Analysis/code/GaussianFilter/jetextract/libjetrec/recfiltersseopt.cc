#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H

#define _XOPEN_SOURCE 600
#include <cstdlib>
#include <jetrec/rec.h>
#include <jetrec/util.h>

namespace jet {

#ifdef HAVE_SSE

	void reconstruction_filtering_t::
	evaluate_step_gaussian_sse(
		float &perp, float gradient[], float hessian[],
		const jet_t &jet,
		const std::vector<track_t>::const_iterator &track_begin,
		const std::vector<track_t>::const_iterator &track_end,
		const collision_geometry_t &geometry) const
	{
#if 0
		float expected_perp;
		float expected_gradient[2];
		float expected_hessian[3];

		evaluate_step_gaussian(
			expected_perp, expected_gradient, expected_hessian, jet,
			track_begin, track_end);
#endif

		/////////////////////////////////////////////////////////////
		// Constants for the azimuthal range reduction
		/////////////////////////////////////////////////////////////
		// M_PI for the SSE2 based azimuthal range reduction
		static float azimuthal_range_reduction_inverse_pi[4]
			__attribute__ ((aligned(16))) = {
			1 / M_PI, 1 / M_PI, 1 / M_PI, 1 / M_PI
		};
		// Integer 1 for the SSE2 based azimuthal range reduction
		static int azimuthal_range_reduction_int_one[4]
			__attribute__ ((aligned(16))) = {
			1, 1, 1, 1
		};
		// -M_PI for both SSE and SSE2 based range reductions. Note
		// that in the SSE version, the sign for the implementation
		// with 1 less register footprint is determined by the
		// possible cmpps ordering, while in SSE2 it does not make any
		// difference (by switching addps <-> subps).
		static float azimuthal_range_reduction_minus_pi[8]
			__attribute__ ((aligned(16))) = {
			-3.14159274F, -3.14159274F,
			-3.14159274F, -3.14159274F,
#ifdef HAVE_SSE2
			8.742278e-8F, 8.742278e-8F,
			8.742278e-8F, 8.742278e-8F
#else // HAVE_SSE2
			// SSE uses compare and subtract by 2 pi, and while the
			// first coefficient is generated arithmetically from the
			// already loaded -pi, the second coefficient is doubled
			// natively for efficiency
			-1.7484556e-7F, -1.7484556e-7F,
			-1.7484556e-7F, -1.7484556e-7F
#endif // HAVE_SSE2
		};
		/////////////////////////////////////////////////////////////
		// Constants for the vectorized exponential function
		/////////////////////////////////////////////////////////////
		// When x >= x_overflow the result is larger than
		// representable in IEEE 754 single precision,
		// vexpf_x_overflow[i] = 128 * M_LN2 for i = 0 .. vl - 1
		static float vexpf_x_overflow[4]
			__attribute__ ((aligned(16))) = {
			88.7228394F, 88.7228394F,
			88.7228394F, 88.7228394F
		};
		// When x < vexpf_x_underflow the result is smaller than
		// representable in IEEE 754 single precision, including
		// denormal, vexpf_x_underflow[i] = -150 M_LN2 for i = 0 .. vl
		// - 1
		static float vexpf_x_underflow[4]
			__attribute__ ((aligned(16))) = {
			-103.972076F, -103.972076F,
			-103.972076F, -103.972076F
		};
		// Range redocution scaling, vexpf_inverse_ln_2[i] = 1 / M_LN2
		// for i = 0 .. vl - 1
		static float vexpf_inverse_ln_2[4]
			__attribute__ ((aligned(16))) = {
			1.44269502F, 1.44269502F,
			1.44269502F, 1.44269502F
		};
		// Software rounding additive bias, vexpf_round_bias[i] = 1.5
		// * pow(2, true_mantissa_length) for i = 0 .. vl - 1
		static float vexpf_round_bias[4]
			__attribute__ ((aligned(16))) = {
			(float)(3 << 22), (float)(3 << 22),
			(float)(3 << 22), (float)(3 << 22)
		};
		// M_LN2 expanded for the Cody and Waite range reduction
		static float vexpf_cody_waite_ln_2[12]
			__attribute__ ((aligned(16))) = {
			0.693145752F, 0.693145752F,
			0.693145752F, 0.693145752F,
			1.42861973e-6F, 1.42861973e-6F,
			1.42861973e-6F, 1.42861973e-6F,
			-1.29053203e-11F, -1.29053203e-11F,
			-1.29053203e-11F, -1.29053203e-11F
		};
		// 7th order minimax polynomial approximation to exp(x) for
		// -M_LN2 <= x <= M_LN2
		static float vexpf_minimax_coefficient[28]
			__attribute__ ((aligned(16))) = {
			2.01509436e-4F, 2.01509436e-4F,
			2.01509436e-4F, 2.01509436e-4F,
			1.41182414e-3F, 1.41182414e-3F,
			1.41182414e-3F, 1.41182414e-3F,
			8.33222549e-3F, 8.33222549e-3F,
			8.33222549e-3F, 8.33222549e-3F,
			4.16603200e-2F, 4.16603200e-2F,
			4.16603200e-2F, 4.16603200e-2F,
			1.66666791e-1F, 1.66666791e-1F,
			1.66666791e-1F, 1.66666791e-1F,
			5.00000477e-1F, 5.00000477e-1F,
			5.00000477e-1F, 5.00000477e-1F,
			1.0F, 1.0F, 1.0F, 1.0F
		};
#ifdef HAVE_SSE2
		// For efficiency reasons, software rounding compares the
		// minimum exponent before removing the rounding bias
		static int vexpf_denormal_exponent[4]
			__attribute__ ((aligned(16))) = {
			0x4b400000 - 125, 0x4b400000 - 125,
			0x4b400000 - 125, 0x4b400000 - 125
		};
		// Denormal exponent offset to compensate for the denormal
		// scaling (below). denormal_exponent_shift[i] =
		// true_mantissa_length for i = 0 .. vl - 1. In the SSE-only
		// implementation, these values are passed as immediates (and
		// are therefore omitted here).
		static int vexpf_denormal_exponent_shift[4]
			__attribute__ ((aligned(16))) = {
			24, 24, 24, 24
		};
#endif // HAVE_SSE2
		// Scaling for denormal results, vexpf_denormal_scale[i] =
		// 1.0F / (float)(1 << vexpf_denormal_exponent_shift[i]) =
		// pow(2, -true_mantissa_length) for i = 0 .. vl - 1.
		static float vexpf_denormal_scale[4]
			__attribute__ ((aligned(16))) = {
			1.0F / (float)(1 << 24), 1.0F / (float)(1 << 24),
			1.0F / (float)(1 << 24), 1.0F / (float)(1 << 24)
		};
		// Maximum float for overflow, which added by itself would
		// give the correctly rounded "infinity" (i.e. infinity for
		// rounding to the nearest (even) and toward +infinity, and
		// 3.40282347e+38 for rounding toward -infinity and 0).
		static float vexpf_float_max[4]
			__attribute__ ((aligned(16))) = {
			3.40282347e+38F, 3.40282347e+38F,
			3.40282347e+38F, 3.40282347e+38F
		};
		// Quiet NaN for masking the result for quiet and signaling
		// NaN arguments
		static unsigned int vexpf_qnan[4]
			__attribute__ ((aligned(16))) = {
			0x7fc00000, 0x7fc00000,
			0x7fc00000, 0x7fc00000
		};
		// SSE requires an additional buffer for %xmm to general
		// purpose register transfer
		float vexpf_buffer[4] __attribute__ ((aligned(16)));
		/////////////////////////////////////////////////////////////

		float scale =
			-0.5F / (_standard_deviation * _standard_deviation);
		float *jet_momentum = const_cast<float *>(jet._momentum._x);
		float deviation_buffer[8] __attribute__ ((aligned(16)));

		__asm__ __volatile__ (
			// Zero %xmm5, %xmm6 for accumulation
			"xorps	%%xmm5, %%xmm5\n\t"
			"xorps	%%xmm6, %%xmm6\n\t"
			// Load the jet momentum into %xmm7[127..64]
			"movss	%2, %%xmm7\n\t"
			// Load the scale into %xmm7[31..0]
			"movhps	%3, %%xmm7"
			: "=m" (scale), "=m" (jet_momentum[2])
			: "m" (scale), "m" (jet_momentum[2])
			: "%xmm5", "%xmm6", "%xmm7");
#ifdef __x86_64__
		__asm__ __volatile__ (
			// Zero %xmm10 .. %xmm15 for accumulation
			"xorps	%%xmm10, %%xmm10\n\t"
			"xorps	%%xmm11, %%xmm11\n\t"
			"xorps	%%xmm12, %%xmm12\n\t"
			"xorps	%%xmm13, %%xmm13\n\t"
			"xorps	%%xmm14, %%xmm14\n\t"
			"xorps	%%xmm15, %%xmm15"
			: :
			: "%xmm10", "%xmm11", "%xmm12", "%xmm13", "%xmm14",
			  "%xmm15");
#endif // __x86_64__

		const long track_size = track_end - track_begin;
		std::vector<track_t>::const_iterator iterator =
			track_begin;

		for(long i = (track_size >> 2); i > 0; i--) {
			IS(__asm__ __volatile__ (
				// Load the jet momentum into %xmm7[127..64]
				"movss	%2, %%xmm7\n\t"
				// Load the scale into %xmm7[31..0]
				"movhps	%3, %%xmm7"
				: "=m" (scale), "=m" (jet_momentum[2])
				: "m" (scale), "m" (jet_momentum[2])
				: "%xmm7"));

			float *track_momentum_0 =
				const_cast<float *>(iterator->_momentum._x);
			iterator++;
			float *track_momentum_1 =
				const_cast<float *>(iterator->_momentum._x);
			iterator++;
			float *track_momentum_2 =
				const_cast<float *>(iterator->_momentum._x);
			iterator++;
			float *track_momentum_3 =
				const_cast<float *>(iterator->_momentum._x);
			iterator++;

			__asm__ __volatile__ (
				"movaps	%%xmm7, %%xmm0\n\t"
				"unpckhps	%%xmm0, %%xmm0\n\t"
				// Load track azimuth and pseudorapidity
				"movlps	%10, %%xmm2\n\t"
				"movlps	%11, %%xmm3\n\t"
				"unpcklps	%%xmm3, %%xmm2\n\t"
				// Now: %xmm0 = [ azimuth[j+1] | azimuth[j] |
				// pseudorapidity[j+1] | pseudorapidity[j] ]
				"movaps	%%xmm0, %%xmm1\n\t"
				"subps	%%xmm2, %%xmm0\n\t"
				//
				"movlps	%12, %%xmm2\n\t"
				"movlps	%13, %%xmm3\n\t"
				"unpcklps	%%xmm3, %%xmm2\n\t"
				// Now: %xmm1 = [ azimuth[j+3] | azimuth[j+2] |
				// pseudorapidity[j+3] | pseudorapidity[j+2] ]
				"subps	%%xmm2, %%xmm1\n\t"
				// %xmm3 is now free
				"movaps	%%xmm1, %%xmm2\n\t"
				"movhlps	%%xmm0, %%xmm1\n\t"
				// Now: %xmm1 = azimuth[j:j+vl-1:1]
				"movlhps	%%xmm2, %%xmm0\n\t"
				// Now: %xmm0 = pseudorapidity[j:j+vl-1:1]
				// %xmm2 is now free
				// Now: %xmm1 = deviation_azimuth[j:j+vl-1:1]
				// %xmm2 is now free
#ifdef HAVE_SSE2
				// SSE2 based azimuthal range reduction
				"movaps	%14, %%xmm2\n\t"
				"mulps	%%xmm1, %%xmm2\n\t"
				"cvttps2dq	%%xmm2, %%xmm2\n\t"
				"movaps	%%xmm2, %%xmm3\n\t"
				// Now: %xmm2 = %xmm3 =
				// (int)(deviation_azimuth[j:j+vl-1:1] / M_PI)
				// Generate bitmask in %xmm3
				"psrad	$32, %%xmm3\n\t"
				// Now: %xmm3 = (%xmm2 >= 0 ? 0 : -1)
				// Load 1 in %xmm4
				// FIXME: Generate using pcmpxx?
				"movaps	%15, %%xmm4\n\t"
				"por	%%xmm4, %%xmm3\n\t"
				// Now: %xmm3 = (%xmm2 >= 0 ? 1 : -1)
				"paddd	%%xmm3, %%xmm2\n\t"
				// Now: %xmm2 = (%xmm2 >= 0 ? %xmm2 + 1 : %xmm2 - 1),
				// i.e. -1 -> -2, 0 -> 1, 1 -> 2
				// %xmm3 is now free
				// Calculate %xmm4 = %xmm2 & ~%xmm4, i.e. for %xmm4 =
				// 1, setting the last bit = 0 and therefore %xmm4 =
				// %xmm2 - %xmm2 % 2 = (%xmm2 / 2) * 2.
				"pandn	%%xmm2, %%xmm4\n\t"
				// Now: -1 -> -2, 0 -> 0, 1 -> 2
				// %xmm2 is now free
				"cvtdq2ps	%%xmm4, %%xmm4\n\t"
				// Add using Cody and Waite range reduction
				"movaps	%16, %%xmm3\n\t"
				"mulps	%%xmm4, %%xmm3\n\t"
				"addps	%%xmm3, %%xmm1\n\t"
				"mulps	%17, %%xmm4\n\t"
				"addps	%%xmm4, %%xmm1\n\t"
#else // HAVE_SSE2
				// SSE based azimuthal range reduction
				"movaps	%%xmm1, %%xmm2\n\t"
				// Load -M_PI
				"movaps	%16, %%xmm3\n\t"
				// Compare deviation_azimuth < -M_PI
				"cmpltps	%%xmm3, %%xmm2\n\t"
				"andps	%%xmm3, %%xmm2\n\t"
				"addps	%%xmm2, %%xmm2\n\t"
				// Now: %xmm2 = deviation_azimuth[j] < -M_PI ? -2 *
				// M_PI : 0, with C_0 = -2 * M_PI for the Cody and
				// Waite range reduction
				"subps	%%xmm2, %%xmm1\n\t"
				// Load C_1 into %xmm4
				"movaps	%17, %%xmm4\n\t"
				"movaps	%%xmm3, %%xmm2\n\t"
				// Apply C_1
				"andps	%%xmm4, %%xmm2\n\t"
				"subps	%%xmm2, %%xmm1\n\t"
				// Generate M_PI in %xmm2
				"xorps	%%xmm2, %%xmm2\n\t"
				"subps	%%xmm3, %%xmm2\n\t"
				// Compare deviation_azimuth > M_PI
				"movaps	%%xmm2, %%xmm3\n\t"
				"cmpltps	%%xmm1, %%xmm3\n\t"
				"addps	%%xmm2, %%xmm2\n\t"
				"andps	%%xmm3, %%xmm2\n\t"
				// Now: %xmm2 = deviation_azimuth[j] > M_PI ? 2 * M_PI
				// : 0
				"subps	%%xmm2, %%xmm1\n\t"
				// Apply C_1, note the sign difference
				"andps	%%xmm4, %%xmm3\n\t"
				"addps	%%xmm3, %%xmm1\n\t"
				// %xmm2, %xmm3 are now free
				// Now: %xmm1 = deviation_azimuth[j:j+vl-1:1]
#endif // HAVE_SSE2
#ifdef __x86_64__
				// AMD64/Intel 64
				"movaps	%%xmm0, %%xmm8\n\t"
				"movaps	%%xmm1, %%xmm9\n\t"
				// Now: %xmm8 = deviation_pseudorapidity[j:j+vl-1:1]
				// Now: %xmm9 = deviation_azimuth[j:j+vl-1:1]
#else // __x86_64__
				// IA-32 is one register short for both
				// deviation_pseudorapidity and deviation_azimuth,
				// therefore store deviation_azimuth in memory.
				"movaps	%%xmm0, %8\n\t"
				"movaps	%%xmm1, %9\n\t"
#endif // __x86_64__
				// Load the exponential scale, moved here due to
				// pipelining
				"movss	%%xmm7, %%xmm2\n\t"
				// Square deviation_pseudorapidity
				"mulps	%%xmm0, %%xmm0\n\t"
				// Now: %xmm0 =
				// deviation_pseudorapidity**2[j:j+vl-1:1]
				"shufps	$0x0, %%xmm2, %%xmm2\n\t"
				// Square deviation_azimuth
				"mulps	%%xmm1, %%xmm1\n\t"
				// Now: %xmm1 = deviation_azimuth**2[j:j+vl-1:1]
				// Accumulate and scale
				"addps	%%xmm1, %%xmm0\n\t"
				"mulps	%%xmm2, %%xmm0"
				// Now: %xmm0 = exp_argument[j:j+vl-1:1]
				: "=m" (track_momentum_0[2]),
				  "=m" (track_momentum_1[2]),
				  "=m" (track_momentum_2[2]),
				  "=m" (track_momentum_3[2]),
				  // 4
				  "=m" (*azimuthal_range_reduction_inverse_pi),
				  "=m" (*azimuthal_range_reduction_int_one),
				  "=m" (azimuthal_range_reduction_minus_pi[0]),
				  "=m" (azimuthal_range_reduction_minus_pi[4]),
				  // 8
				  "=m" (deviation_buffer[0]),
				  "=m" (deviation_buffer[4])
				: // 10
				  "m" (track_momentum_0[2]),
				  "m" (track_momentum_1[2]),
				  "m" (track_momentum_2[2]),
				  "m" (track_momentum_3[2]),
				  // 14
				  "m" (*azimuthal_range_reduction_inverse_pi),
				  "m" (*azimuthal_range_reduction_int_one),
				  "m" (azimuthal_range_reduction_minus_pi[0]),
				  "m" (azimuthal_range_reduction_minus_pi[4]),
				  // 18
				  "m" (deviation_buffer[0]),
				  "m" (deviation_buffer[4])
				: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4", "%xmm7"
#ifdef __x86_64__
				  , "%xmm8", "%xmm9"
#endif // __x86_64__
				);
#ifdef DEBUG_ASSEMBLY
			ID(float debug_buffer[12] __attribute__ ((aligned(16))));
			IS(__asm__ __volatile__ (
				"movaps	%%xmm0, %0\n\t"
				"movaps	%%xmm1, %1\n\t"
				"movaps	%%xmm4, %2"
				: "=m" (debug_buffer[0]), "=m" (debug_buffer[4]),
				  "=m" (debug_buffer[8])
				: "m" (debug_buffer[0]), "m" (debug_buffer[4]),
				  "m" (debug_buffer[8])
				: "%xmm0", "%xmm1", "%xmm2"));
			for(int j = 0; j < 4; j++) {
				ID(const float deviation_pseudorapidity =
				   iterator[j - 4].momentum().pseudorapidity() -
				   jet.momentum().pseudorapidity());
				ID(const float deviation_azimuth =
				   angular_range_reduce(
					iterator[j - 4].momentum().azimuth() -
					jet.momentum().azimuth()));
				ID(const float radial_deviation =
				   deviation_pseudorapidity *
				   deviation_pseudorapidity +
				   deviation_azimuth * deviation_azimuth);
				ID(const float exponential_factor =
				   -radial_deviation /
				   (2 * _standard_deviation * _standard_deviation));
				I(FEQ(debug_buffer[j], exponential_factor));
			}
#endif // DEBUG_ASSEMBLY
			__asm__ __volatile__ (
				// Load vexpf_round_bias[0 : vl - 1 : 1]
				"movaps	%3, %%xmm1\n\t"
				// We need x later for overflow/underflow test
				"movaps	%%xmm0, %%xmm2\n\t"
				// Now: %xmm2 = x[i : i + vl - 1 : 1]
				// Scale by 1 / ln(2)
				"movaps	%%xmm0, %%xmm3\n\t"
				"mulps	%2, %%xmm3\n\t"
				// Step 1 of rounding by bit manipulation: (floating
				// point) add 1.5 * pow(2, true_mantissa_length)
				"addps	%%xmm1, %%xmm3\n\t"
				"movaps	%%xmm3, %%xmm4\n\t"
				// Now: %xmm4 = k[i : i + vl - 1 : 1] + bias
				"subps	%%xmm1, %%xmm3"
				// Now: %xmm3 = fl(k[i : i + vl - 1 : 1])
				: "=m" (*vexpf_inverse_ln_2), "=m" (*vexpf_round_bias)
				: "m" (*vexpf_inverse_ln_2), "m" (*vexpf_round_bias)
				: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4");
			__asm__ __volatile__ (
				// Cody and Waite range reduction with 3 coefficients
				"movaps	%%xmm3, %%xmm1\n\t"
				"mulps	%3, %%xmm1\n\t"
				"subps	%%xmm1, %%xmm0\n\t"
				"movaps	%%xmm3, %%xmm1\n\t"
				"mulps	%4, %%xmm1\n\t"
				"subps	%%xmm1, %%xmm0\n\t"
				"movaps	%%xmm3, %%xmm1\n\t"
				"mulps	%5, %%xmm1\n\t"
				"subps	%%xmm1, %%xmm0"
				// Now: %xmm0 = xs[i : i + vl - 1 : 1]
				: // 0
				  "=m" (vexpf_cody_waite_ln_2[0]),
				  "=m" (vexpf_cody_waite_ln_2[4]),
				  "=m" (vexpf_cody_waite_ln_2[8])
				: // 3
				  "m" (vexpf_cody_waite_ln_2[0]),
				  "m" (vexpf_cody_waite_ln_2[4]),
				  "m" (vexpf_cody_waite_ln_2[8])
				: "%xmm0", "%xmm1", "%xmm3");
			__asm__ __volatile__ (
				// Horner rule minimax polynomial evaluation
				"movaps	%%xmm0, %%xmm1\n\t"
				"mulps	%7, %%xmm1\n\t"
				"addps	%8, %%xmm1\n\t"
				"mulps	%%xmm0, %%xmm1\n\t"
				"addps	%9, %%xmm1\n\t"
				"mulps	%%xmm0, %%xmm1\n\t"
				"addps	%10, %%xmm1\n\t"
				"mulps	%%xmm0, %%xmm1\n\t"
				"addps	%11, %%xmm1\n\t"
				"mulps	%%xmm0, %%xmm1\n\t"
				"addps	%12, %%xmm1\n\t"
				"mulps	%%xmm0, %%xmm1\n\t"
				"mulps	%%xmm0, %%xmm1\n\t"
				"addps	%%xmm1, %%xmm0\n\t"
				"addps	%13, %%xmm0"
				// Now: %xmm0 = g[i : i + vl - 1 : 1]
				: // 0
				  "=m" (vexpf_minimax_coefficient[0]),
				  "=m" (vexpf_minimax_coefficient[4]),
				  "=m" (vexpf_minimax_coefficient[8]),
				  "=m" (vexpf_minimax_coefficient[12]),
				  "=m" (vexpf_minimax_coefficient[16]),
				  "=m" (vexpf_minimax_coefficient[20]),
				  "=m" (vexpf_minimax_coefficient[24])
				: // 7
				  "m" (vexpf_minimax_coefficient[0]),
				  "m" (vexpf_minimax_coefficient[4]),
				  "m" (vexpf_minimax_coefficient[8]),
				  "m" (vexpf_minimax_coefficient[12]),
				  "m" (vexpf_minimax_coefficient[16]),
				  "m" (vexpf_minimax_coefficient[20]),
				  "m" (vexpf_minimax_coefficient[24])
				: "%xmm0", "%xmm1");
#ifdef HAVE_SSE2
			__asm__ __volatile__ (
				// Now: %xmm0 = g[i : i + vl - 1 : 1]
				// Now: %xmm2 = x[i : i + vl - 1 : 1]
				// Now: %xmm4 = k[i : i + vl - 1 : 1] + bias
				// Mask k > vexpf_denormal_exponent
				"movdqa	%3, %%xmm1\n\t"
				"pcmpgtd	%%xmm4, %%xmm1\n\t"
				// Now: %xmm1 = k[i : i + vl - 1 : 1] <
				// vexpf_denormal_exponent ? 0xffffffff : 0
				"movmskps	%%xmm1, %%eax\n\t"
				"testb	$0xf, %%al\n\t"
				"jz	0f\n\t"
				// Load and mask denormal exponent shift
				"movdqa	%4, %%xmm3\n\t"
				"pand	%%xmm1, %%xmm3\n\t"
				// Add denormal exponent shift
				"paddd	%%xmm4, %%xmm3\n\t"
				// %xmm4 is now free
				// Shift exponent to the IEEE 754 single precision
				// exponent position
				"pslld	$23, %%xmm3\n\t"
				// Add exponent
				"paddd	%%xmm3, %%xmm0\n\t"
				// %xmm3 is now free
				// Multiply by denormal scale and mask
				"movaps	%%xmm0, %%xmm3\n\t"
				"mulps	%5, %%xmm3\n\t"
				"andps	%%xmm1, %%xmm3\n\t"
				// Inversely mask the unscaled result
				"andnps	%%xmm0, %%xmm1\n\t"
				// Merge masked results
				"orps	%%xmm3, %%xmm1\n\t"
				// Store
				"movaps	%%xmm1, %%xmm0\n\t"
				"jmp	1f\n"
				"0:\n\t"
				// Same as above, without denormal handling
				"pslld	$23, %%xmm4\n\t"
				"paddd	%%xmm4, %%xmm0\n"
				"1:\n\t"
				// Now: %xmm0 = exp_x[j : j + vl - 1 : 1]
				: "=m" (*vexpf_denormal_exponent),
				  "=m" (*vexpf_denormal_exponent_shift),
				  "=m" (*vexpf_denormal_scale)
				: "m" (*vexpf_denormal_exponent),
				  "m" (*vexpf_denormal_exponent_shift),
				  "m" (*vexpf_denormal_scale)
				: "%eax", "%xmm0", "%xmm1", "%xmm3", "%xmm4");
#else // HAVE_SSE2
			// Apply exponents in %eax, with denormal handling
			__asm__ __volatile__ (
				// Now: %xmm0 = g[i : i + vl - 1 : 1]
				// Now: %xmm2 = x[i : i + vl - 1 : 1]
				// Now: %xmm4 = k[i : i + vl - 1 : 1] + bias
				// Store g[i : i + vl - 1 : 1] into the output, using
				// it as temprary storage
				"movaps	%%xmm0, %0\n\t"
				// Load k[0] into %eax
				"pextrw	$0x0, %%xmm4, %%eax\n\t"
				// Test if result is denormal
				"cmpw	$-126, %%ax\n\t"
				"jg	0f\n\t"
				// Result is denormal:
				// Prescale exponent by the effective mantissa length
				"addw	$24, %%ax\n\t"
				// Shift exponent to the floating point exponent
				// position
				"shll	$23, %%eax\n\t"
				// Scale g by the exponent
				"addl	%%eax, %0\n\t"
				// Rescale to obtain the denormal result
				"movss	%2, %%xmm0\n\t"
				"mulss	%3, %%xmm0\n\t"
				"movss	%%xmm0, %0\n"
				"0:\n\t"
				// Result is normal:
				// Analogously to the denormal case, without scaling
				"shll	$23, %%eax\n\t"
				"addl	%%eax, %0"
				: "=m" (vexpf_buffer[0]),
				  "=m" (*vexpf_denormal_scale)
				: "m" (vexpf_buffer[0]),
				  "m" (*vexpf_denormal_scale)
				: "%eax", "%xmm0", "%xmm4");
			// Repeat for k[i + 1 : i + vl - 1 : 1]
			__asm__ __volatile__ (
				"pextrw	$0x2, %%xmm4, %%eax\n\t"
				"cmpw	$-126, %%ax\n\t"
				"jg	0f\n\t"
				"addw	$24, %%ax\n\t"
				"shll	$23, %%eax\n\t"
				"addl	%%eax, %0\n\t"
				"movss	%2, %%xmm0\n\t"
				"mulss	%3, %%xmm0\n\t"
				"movss	%%xmm0, %0\n"
				"0:\n\t"
				"shll	$23, %%eax\n\t"
				"addl	%%eax, %0"
				: "=m" (vexpf_buffer[1]),
				  "=m" (*vexpf_denormal_scale)
				: "m" (vexpf_buffer[1]),
				  "m" (*vexpf_denormal_scale)
				: "%eax", "%xmm0", "%xmm4");
			__asm__ __volatile__ (
				"pextrw	$0x4, %%xmm4, %%eax\n\t"
				"cmpw	$-126, %%ax\n\t"
				"jg	0f\n\t"
				"addw	$24, %%ax\n\t"
				"shll	$23, %%eax\n\t"
				"addl	%%eax, %0\n\t"
				"movss	%2, %%xmm0\n\t"
				"mulss	%3, %%xmm0\n\t"
				"movss	%%xmm0, %0\n"
				"0:\n\t"
				"shll	$23, %%eax\n\t"
				"addl	%%eax, %0"
				: "=m" (vexpf_buffer[2]),
				  "=m" (*vexpf_denormal_scale)
				: "m" (vexpf_buffer[2]),
				  "m" (*vexpf_denormal_scale)
				: "%eax", "%xmm0", "%xmm4");
			__asm__ __volatile__ (
				"pextrw	$0x6, %%xmm4, %%eax\n\t"
				"cmpw	$-126, %%ax\n\t"
				"jg	0f\n\t"
				"addw	$24, %%ax\n\t"
				"shll	$23, %%eax\n\t"
				"addl	%%eax, %0\n\t"
				"movss	%2, %%xmm0\n\t"
				"mulss	%3, %%xmm0\n\t"
				"movss	%%xmm0, %0\n"
				"0:\n\t"
				"shll	$23, %%eax\n\t"
				"addl	%%eax, %0"
				: "=m" (vexpf_buffer[3]),
				  "=m" (*vexpf_denormal_scale)
				: "m" (vexpf_buffer[3]),
				  "m" (*vexpf_denormal_scale)
				: "%eax", "%xmm0", "%xmm4");
#endif // HAVE_SSE2
			__asm__ __volatile__ (
				// Now: %xmm0 = exp_x[i : i + vl - 1 : 1]
				// Now: %xmm2 = x[i : i + vl - 1 : 1]
				"movaps	%%xmm2, %%xmm3\n\t"
				// Test x < vexpf_x_underflow
				"cmpltps	%5, %%xmm3\n\t"
				// Zero underflows by simple masking
#ifdef HAVE_SSE2
				"andnps	%%xmm0, %%xmm3\n\t"
#else // HAVE_SSE2
				"andnps	%7, %%xmm3\n\t"
#endif // HAVE_SSE2
				// Now: %xmm3 = exp_x[i : i + vl - 1 : 1]
				// Test x >= vexpf_x_overflow
				"movaps	%6, %%xmm1\n\t"
				"cmpleps	%%xmm2, %%xmm1\n\t"
				// Generating rounded "infinity" is expensive, and
				// jumping over for vector nonoverflow becomes
				// efficient.
				"movmskps	%%xmm1, %%eax\n\t"
				"testb	$0xf, %%al\n\t"
				"jz	0f\n\t"
				// Set overflows to the rounded "infinity"
				"movaps	%%xmm1, %%xmm4\n\t"
				// Generate rounded "infinity"
				"movaps	%8, %%xmm0\n\t"
				"addps	%%xmm0, %%xmm0\n\t"
				// Merge by masking
				"andps	%%xmm0, %%xmm4\n\t"
				"andnps	%%xmm3, %%xmm1\n\t"
				"orps	%%xmm1, %%xmm4\n\t"
				"movaps	%%xmm4, %%xmm3\n"
				"0:\n\t"
				// Test x == (Q, S)NaN
				"movaps	%%xmm2, %%xmm0\n\t"
				"cmpeqps	%%xmm0, %%xmm0\n\t"
				"andps	%%xmm0, %%xmm3\n\t"
				"orps	%9, %%xmm2\n\t"
				"andnps	%%xmm2, %%xmm0\n\t"
				"orps	%%xmm3, %%xmm0"
				// Now: %xmm0 = exp_x[i : i + vl - 1 : 1]
				: "=m" (*vexpf_x_underflow),
				  "=m" (*vexpf_x_overflow),
				  "=m" (*vexpf_buffer), "=m" (*vexpf_float_max),
				  "=m" (*vexpf_qnan)
				: // 5
				  "m" (*vexpf_x_underflow),
				  "m" (*vexpf_x_overflow),
				  "m" (*vexpf_buffer), "m" (*vexpf_float_max),
				  "m" (*vexpf_qnan)
				: "%eax", "%xmm0", "%xmm1", "%xmm2", "%xmm3",
				  "%xmm4");
			// Generate perp[j:j+vl-1:1], gradient[j:j+vl-1:1][0:1],
			// hessian_sum[j:j+vl-1:1][1]
			__asm__ __volatile__ (
				// Load the transverse momenta
				"movss	%6, %%xmm1\n\t"
				"movss	%7, %%xmm2\n\t"
				"unpcklps	%%xmm2, %%xmm1\n\t"
				"movss	%8, %%xmm2\n\t"
				"movss	%9, %%xmm3\n\t"
				"unpcklps	%%xmm3, %%xmm2\n\t"
				// %xmm3 is now free
				"movlhps	%%xmm2, %%xmm1\n\t"
				// %xmm2 is now free
				// Apply the transverse momenta to the result from
				// exponentiation
				"mulps	%%xmm1, %%xmm0\n\t"
				// Now: %xmm0 = perp[j:j+vl-1:1]
				// %xmm1 is now free
				// Load the exponential scale
				"movss	%%xmm7, %%xmm1\n\t"
				"addss	%%xmm1, %%xmm1\n\t"
				"shufps	$0x0, %%xmm1, %%xmm1\n\t"
				"movaps	%%xmm1, %%xmm3\n\t"
				"mulps	%%xmm3, %%xmm3\n\t"
				"mulps	%%xmm0, %%xmm1\n\t"
				"movaps	%%xmm1, %%xmm2\n\t"
				"mulps	%%xmm0, %%xmm3\n\t"
#ifdef __x86_64__
				// AMD64/Intel 64
				// Now: %xmm8 = deviation_pseudorapidity[j:j+vl-1:1]
				// Now: %xmm9 = deviation_azimuth[j:j+vl-1:1]
				"addps	%%xmm0, %%xmm10\n\t"
				"mulps	%%xmm8, %%xmm1\n\t"
				"addps	%%xmm1, %%xmm11\n\t"
				"mulps	%%xmm9, %%xmm2\n\t"
				"addps	%%xmm2, %%xmm12\n\t"
				"mulps	%%xmm8, %%xmm3\n\t"
				"mulps	%%xmm9, %%xmm3\n\t"
				"addps	%%xmm3, %%xmm13\n\t"
#else // __x86_64__
				// IA-32
				"mulps	%10, %%xmm1\n\t"
				"mulps	%11, %%xmm2\n\t"
				"mulps	%10, %%xmm3\n\t"
				"mulps	%11, %%xmm3"
#endif // __x86_64__
				// Hessians are stored as {h(1, 1), h(2, 2), h(1, 2)}
				// Now: %xmm0 = perp[j:j+vl-1:1][0]
				// Now: %xmm1 = gradient[j:j+vl-1:1][0]
				// Now: %xmm2 = gradient[j:j+vl-1:1][1]
				// Now: %xmm3 = hessian[j:j+vl-1:1][1]
				: "=m" (track_momentum_0[1]),
				  "=m" (track_momentum_1[1]),
				  "=m" (track_momentum_2[1]),
				  "=m" (track_momentum_3[1]),
				  "=m" (deviation_buffer[0]),
				  "=m" (deviation_buffer[4])
				: "m" (track_momentum_0[1]),
				  "m" (track_momentum_1[1]),
				  "m" (track_momentum_2[1]),
				  "m" (track_momentum_3[1]),
				  "m" (deviation_buffer[0]),
				  "m" (deviation_buffer[4])
				: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
				  "%xmm7"
#ifdef __x86_64__
				  , "%xmm8", "%xmm9", "%xmm10", "%xmm11", "%xmm12",
				  "%xmm13", "%xmm14", "%xmm15"
#endif // __x86_64__
				);
#ifdef DEBUG_ASSEMBLY
			ID(float debug_buffer[12] __attribute__ ((aligned(16))));
			IS(__asm__ __volatile__ (
				"movaps	%%xmm0, %0\n\t"
				"movaps	%%xmm1, %1\n\t"
				"movaps	%%xmm2, %2"
				: "=m" (debug_buffer[0]), "=m" (debug_buffer[4]),
				  "=m" (debug_buffer[8])
				: "m" (debug_buffer[0]), "m" (debug_buffer[4]),
				  "m" (debug_buffer[8])
				: "%xmm0", "%xmm1", "%xmm2"));
			for(int j = 0; j < 4; j++) {
				ID(const float deviation_pseudorapidity =
				   iterator[j - 4].momentum().pseudorapidity() -
				   jet.momentum().pseudorapidity());
				ID(const float deviation_azimuth =
				   angular_range_reduce(
					iterator[j - 4].momentum().azimuth() -
					jet.momentum().azimuth()));
				ID(const float radial_deviation =
				   deviation_pseudorapidity *
				   deviation_pseudorapidity +
				   deviation_azimuth * deviation_azimuth);
				ID(const float exponential_factor =
				   iterator[j - 4].momentum().perp() *
				   exp(-radial_deviation /
					   (2 * _standard_deviation *
						_standard_deviation)));
				I(FEQ(debug_buffer[j], exponential_factor));
			}
#endif // DEBUG_ASSEMBLY
#ifndef __x86_64__
			__asm__ __volatile__ (
#ifdef HAVE_SSE3
				// SSE3 (in place) horizontal add
				"haddps	%%xmm1, %%xmm1\n\t"
				"haddps	%%xmm1, %%xmm1\n\t"
				// %xmm1 -> %xmm2
				"haddps	%%xmm2, %%xmm2\n\t"
				"haddps	%%xmm2, %%xmm2\n\t"
				// %xmm2 -> %xmm3
				"haddps	%%xmm3, %%xmm3\n\t"
				"haddps	%%xmm3, %%xmm3\n\t"
#else // HAVE_SSE3
				// SSE (in place, simulated) horizontal add
				"movhlps	%%xmm1, %%xmm4\n\t"
				"addps	%%xmm1, %%xmm4\n\t"
				"movaps	%%xmm4, %%xmm1\n\t"
				"shufps	$0x55, %%xmm1, %%xmm1\n\t"
				"addss	%%xmm4, %%xmm1\n\t"
				// %xmm1 -> %xmm2
				"movhlps	%%xmm2, %%xmm4\n\t"
				"addps	%%xmm2, %%xmm4\n\t"
				"movaps	%%xmm4, %%xmm2\n\t"
				"shufps	$0x55, %%xmm2, %%xmm2\n\t"
				"addss	%%xmm4, %%xmm2\n\t"
				// %xmm2 -> %xmm3
				"movhlps	%%xmm3, %%xmm4\n\t"
				"addps	%%xmm3, %%xmm4\n\t"
				"movaps	%%xmm4, %%xmm3\n\t"
				"shufps	$0x55, %%xmm3, %%xmm3\n\t"
				"addss	%%xmm4, %%xmm3\n\t"
#endif // HAVE_SSE3
				// Pack %xmm1, %xmm2
				"unpcklps	%%xmm2, %%xmm1\n\t"
				"movlhps	%%xmm1, %%xmm4\n\t"
				// %xmm1, %xmm2 are now free
#ifdef HAVE_SSE3
				// SSE3 (out of place) horizontal add
				"movaps	%%xmm0, %%xmm1\n\t"
				"haddps	%%xmm1, %%xmm1\n\t"
				"haddps	%%xmm1, %%xmm1\n\t"
#else // HAVE_SSE3
				// SSE (out of place) horizontal add
				"movhlps	%%xmm0, %%xmm1\n\t"
				"addps	%%xmm0, %%xmm1\n\t"
				"movaps	%%xmm1, %%xmm2\n\t"
				"shufps	$0x55, %%xmm2, %%xmm2\n\t"
				"addss	%%xmm2, %%xmm1\n\t"
#endif // HAVE_SSE3
				"movss	%%xmm1, %%xmm4\n\t"
				// Now: %xmm4 = [ gradient_sum[1] | gradient_sum[0] |
				// * | perp_sum ]
				"addps	%%xmm4, %%xmm5"
				: : : "%xmm0", "%xmm1", "%xmm2", "%xmm4", "%xmm5");
#endif // __x86_64__

			static float one[4] __attribute__ ((aligned(16))) = {
				1.0F, 1.0F, 1.0F, 1.0F
			};

			// Generate hessian_sum[0:2:2]
#ifdef __x86_64__
			__asm__ __volatile__ (
				// AMD64/Intel 64
				// Now: %xmm8 = deviation_pseudorapidity[j:j+vl-1:1]
				// Now: %xmm9 = deviation_azimuth[j:j+vl-1:1]
				// Load the exponential scale
				"movss	%%xmm7, %%xmm1\n\t"
				"addss	%%xmm1, %%xmm1\n\t"
				"shufps	$0x0, %%xmm1, %%xmm1\n\t"
				// Now: %xmm1 = (2 * scale)[j:j+vl-1:1]
				"mulps	%%xmm8, %%xmm8\n\t"
				"mulps	%%xmm1, %%xmm8\n\t"
				"addps	%1, %%xmm8\n\t"
				"mulps	%%xmm1, %%xmm8\n\t"
				"mulps	%%xmm0, %%xmm8\n\t"
				"addps	%%xmm8, %%xmm14\n\t"
				// %xmm8 -> %xmm9
				"mulps	%%xmm9, %%xmm9\n\t"
				"mulps	%%xmm1, %%xmm9\n\t"
				"addps	%1, %%xmm9\n\t"
				"mulps	%%xmm1, %%xmm9\n\t"
				"mulps	%%xmm0, %%xmm9\n\t"
				"addps	%%xmm9, %%xmm15\n\t"
				: "=m" (*one)
				: "m" (*one)
				: "%xmm0", "%xmm8", "%xmm9", "%xmm14", "%xmm15");
#else // __x86_64__
			__asm__ __volatile__ (
				// Load the exponential scale
				"movss	%%xmm7, %%xmm1\n\t"
				"addss	%%xmm1, %%xmm1\n\t"
				"shufps	$0x0, %%xmm1, %%xmm1\n\t"
#ifdef __x86_64__
				// AMD64/Intel 64
				// Now: %xmm8 = deviation_pseudorapidity[j:j+vl-1:1]
				// Now: %xmm9 = deviation_azimuth[j:j+vl-1:1]
				"mulps	%%xmm8, %%xmm8\n\t"
				"mulps	%%xmm1, %%xmm8\n\t"
				"addps	%3, %%xmm8\n\t"
				"mulps	%%xmm1, %%xmm8\n\t"
				"mulps	%%xmm0, %%xmm8\n\t"
				// %xmm8 -> %xmm9
				"mulps	%%xmm9, %%xmm9\n\t"
				"mulps	%%xmm1, %%xmm9\n\t"
				"addps	%3, %%xmm9\n\t"
				"mulps	%%xmm1, %%xmm9\n\t"
				"mulps	%%xmm0, %%xmm9\n\t"
				// %xmm1 is now free
				// Now: %xmm6 = [ * | hessian_sum[2] |
				// hessian_sum[1] | hessian_sum[0] ]
#ifdef HAVE_SSE3
				// SSE3 (in place) horizontal add
				"haddps	%%xmm8, %%xmm8\n\t"
				"haddps	%%xmm8, %%xmm8\n\t"
				// %xmm8 -> %xmm9
				"haddps	%%xmm9, %%xmm9\n\t"
				"haddps	%%xmm9, %%xmm9\n\t"
#else // HAVE_SSE3
				// SSE (in place, simulated) horizontal add
				"movhlps	%%xmm8, %%xmm1\n\t"
				"addps	%%xmm8, %%xmm1\n\t"
				"movaps	%%xmm1, %%xmm8\n\t"
				"shufps	$0x55, %%xmm8, %%xmm8\n\t"
				"addss	%%xmm1, %%xmm8\n\t"
				// %xmm8 -> %xmm9
				"movhlps	%%xmm9, %%xmm1\n\t"
				"addps	%%xmm9, %%xmm1\n\t"
				"movaps	%%xmm1, %%xmm9\n\t"
				"shufps	$0x55, %%xmm9, %%xmm9\n\t"
				"addss	%%xmm1, %%xmm9\n\t"
#endif // HAVE_SSE3
				"unpcklps	%%xmm9, %%xmm8\n\t"
				"movlhps	%%xmm3, %%xmm8\n\t"
				"addps	%%xmm8, %%xmm6"
#else // __x86_64__
				// IA-32
				"movaps	%4, %%xmm2\n\t"
				// Now same as AMD64/Intel 64, with %xmm8 -> %xmm2
				"mulps	%%xmm2, %%xmm2\n\t"
				"mulps	%%xmm1, %%xmm2\n\t"
				"addps	%3, %%xmm2\n\t"
				"mulps	%%xmm1, %%xmm2\n\t"
				"mulps	%%xmm0, %%xmm2\n\t"
				// %xmm2 -> %xmm4
				"movaps	%5, %%xmm4\n\t"
				"mulps	%%xmm4, %%xmm4\n\t"
				"mulps	%%xmm1, %%xmm4\n\t"
				"addps	%3, %%xmm4\n\t"
				"mulps	%%xmm1, %%xmm4\n\t"
				"mulps	%%xmm0, %%xmm4\n\t"
#ifdef HAVE_SSE3
				// SSE3 (in place) horizontal add
				"haddps	%%xmm2, %%xmm2\n\t"
				"haddps	%%xmm2, %%xmm2\n\t"
				// %xmm8 -> %xmm9
				"haddps	%%xmm4, %%xmm4\n\t"
				"haddps	%%xmm4, %%xmm4\n\t"
#else // HAVE_SSE3
				// SSE (in place, simulated) horizontal add
				"movhlps	%%xmm2, %%xmm1\n\t"
				"addps	%%xmm2, %%xmm1\n\t"
				"movaps	%%xmm1, %%xmm2\n\t"
				"shufps	$0x55, %%xmm2, %%xmm2\n\t"
				"addss	%%xmm1, %%xmm2\n\t"
				// %xmm2 -> %xmm4
				"movhlps	%%xmm4, %%xmm1\n\t"
				"addps	%%xmm4, %%xmm1\n\t"
				"movaps	%%xmm1, %%xmm4\n\t"
				"shufps	$0x55, %%xmm4, %%xmm4\n\t"
				"addss	%%xmm1, %%xmm4\n\t"
#endif // HAVE_SSE3
				"unpcklps	%%xmm4, %%xmm2\n\t"
				"movlhps	%%xmm3, %%xmm2\n\t"
				"addps	%%xmm2, %%xmm6"
#endif // __x86_64__
				: "=m" (*one), "=m" (deviation_buffer[0]),
				  "=m" (deviation_buffer[4])
				: "m" (*one), "m" (deviation_buffer[0]),
				  "m" (deviation_buffer[4])
				: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
				  "%xmm6", "%xmm7"
#ifdef __x86_64__
				  , "%xmm8", "%xmm9"
#endif // __x86_64__
				);
#endif // __x86_64__
		}

		const unsigned int track_size_mod_4 = track_size & 3;

		if(track_size_mod_4 != 0) {
#ifdef __INTEL_COMPILER
			float *track_momentum_0;
			float *track_momentum_1;
			float *track_momentum_2;
#else // __INTEL_COMPILER
			// GCC is extremely stupid when analyzing what is
			// initialized...
			float *track_momentum_0 = NULL;
			float *track_momentum_1 = NULL;
			float *track_momentum_2 = NULL;
#endif // __INTEL_COMPILER
			static float track_momentum_3[4]
				__attribute__ ((aligned(16))) = {
				0.0F, 0.0F, 0.0F, 0.0F
			};

			switch(track_size_mod_4) {
			case 1:
				track_momentum_0 =
					const_cast<float *>(iterator->_momentum._x);
				track_momentum_1 = track_momentum_3;	// zero
				track_momentum_2 = track_momentum_3;	// zero
				break;
			case 2:
				track_momentum_0 =
					const_cast<float *>(iterator++->_momentum._x);
				track_momentum_1 =
					const_cast<float *>(iterator->_momentum._x);
				track_momentum_2 = track_momentum_3;	// zero
				break;
			case 3:
				track_momentum_0 =
					const_cast<float *>(iterator++->_momentum._x);
				track_momentum_1 =
					const_cast<float *>(iterator++->_momentum._x);
				track_momentum_2 =
					const_cast<float *>(iterator->_momentum._x);
				break;
			}

			__asm__ __volatile__ (
				"movaps	%%xmm7, %%xmm0\n\t"
				"unpckhps	%%xmm0, %%xmm0\n\t"
				"movlps	%10, %%xmm2\n\t"
				"movlps	%11, %%xmm3\n\t"
				"unpcklps	%%xmm3, %%xmm2\n\t"
				"movaps	%%xmm0, %%xmm1\n\t"
				"subps	%%xmm2, %%xmm0\n\t"
				"movlps	%12, %%xmm2\n\t"
				"movlps	%13, %%xmm3\n\t"
				"unpcklps	%%xmm3, %%xmm2\n\t"
				"subps	%%xmm2, %%xmm1\n\t"
				"movaps	%%xmm1, %%xmm2\n\t"
				"movhlps	%%xmm0, %%xmm1\n\t"
				"movlhps	%%xmm2, %%xmm0\n\t"
#ifdef HAVE_SSE2
				"movaps	%14, %%xmm2\n\t"
				"mulps	%%xmm1, %%xmm2\n\t"
				"cvttps2dq	%%xmm2, %%xmm2\n\t"
				"movaps	%%xmm2, %%xmm3\n\t"
				"psrad	$32, %%xmm3\n\t"
				"movaps	%15, %%xmm4\n\t"
				"por	%%xmm4, %%xmm3\n\t"
				"paddd	%%xmm3, %%xmm2\n\t"
				"pandn	%%xmm2, %%xmm4\n\t"
				"cvtdq2ps	%%xmm4, %%xmm4\n\t"
				"movaps	%16, %%xmm3\n\t"
				"mulps	%%xmm4, %%xmm3\n\t"
				"addps	%%xmm3, %%xmm1\n\t"
				"mulps	%17, %%xmm4\n\t"
				"addps	%%xmm4, %%xmm1\n\t"
#else // HAVE_SSE2
				"movaps	%%xmm1, %%xmm2\n\t"
				"movaps	%16, %%xmm3\n\t"
				"cmpltps	%%xmm3, %%xmm2\n\t"
				"andps	%%xmm3, %%xmm2\n\t"
				"addps	%%xmm2, %%xmm2\n\t"
				"subps	%%xmm2, %%xmm1\n\t"
				"movaps	%17, %%xmm4\n\t"
				"movaps	%%xmm3, %%xmm2\n\t"
				"andps	%%xmm4, %%xmm2\n\t"
				"subps	%%xmm2, %%xmm1\n\t"
				"xorps	%%xmm2, %%xmm2\n\t"
				"subps	%%xmm3, %%xmm2\n\t"
				"movaps	%%xmm2, %%xmm3\n\t"
				"cmpltps	%%xmm1, %%xmm3\n\t"
				"addps	%%xmm2, %%xmm2\n\t"
				"andps	%%xmm3, %%xmm2\n\t"
				"subps	%%xmm2, %%xmm1\n\t"
				"andps	%%xmm4, %%xmm3\n\t"
				"addps	%%xmm3, %%xmm1\n\t"
#endif // HAVE_SSE2
#ifdef __x86_64__
				"movaps	%%xmm0, %%xmm8\n\t"
				"movaps	%%xmm1, %%xmm9\n\t"
#else // __x86_64__
				"movaps	%%xmm0, %8\n\t"
				"movaps	%%xmm1, %9\n\t"
#endif // __x86_64__
				"movss	%%xmm7, %%xmm2\n\t"
				"mulps	%%xmm0, %%xmm0\n\t"
				"shufps	$0x0, %%xmm2, %%xmm2\n\t"
				"mulps	%%xmm1, %%xmm1\n\t"
				"addps	%%xmm1, %%xmm0\n\t"
				"mulps	%%xmm2, %%xmm0"
				: "=m" (track_momentum_0[2]),
				  "=m" (track_momentum_1[2]),
				  "=m" (track_momentum_2[2]),
				  "=m" (track_momentum_3[2]),
				  "=m" (*azimuthal_range_reduction_inverse_pi),
				  "=m" (*azimuthal_range_reduction_int_one),
				  "=m" (azimuthal_range_reduction_minus_pi[0]),
				  "=m" (azimuthal_range_reduction_minus_pi[4]),
				  "=m" (deviation_buffer[0]),
				  "=m" (deviation_buffer[4])
				: "m" (track_momentum_0[2]),
				  "m" (track_momentum_1[2]),
				  "m" (track_momentum_2[2]),
				  "m" (track_momentum_3[2]),
				  "m" (*azimuthal_range_reduction_inverse_pi),
				  "m" (*azimuthal_range_reduction_int_one),
				  "m" (azimuthal_range_reduction_minus_pi[0]),
				  "m" (azimuthal_range_reduction_minus_pi[4]),
				  "m" (deviation_buffer[0]),
				  "m" (deviation_buffer[4])
				: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm7"
#ifdef __x86_64__
				  , "%xmm8", "%xmm9"
#endif // __x86_64__
				);
			__asm__ __volatile__ (
				"movaps	%3, %%xmm1\n\t"
				"movaps	%%xmm0, %%xmm2\n\t"
				"movaps	%%xmm0, %%xmm3\n\t"
				"mulps	%2, %%xmm3\n\t"
				"addps	%%xmm1, %%xmm3\n\t"
				"movaps	%%xmm3, %%xmm4\n\t"
				"subps	%%xmm1, %%xmm3"
				: "=m" (*vexpf_inverse_ln_2), "=m" (*vexpf_round_bias)
				: "m" (*vexpf_inverse_ln_2), "m" (*vexpf_round_bias)
				: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4");
			__asm__ __volatile__ (
				"movaps	%%xmm3, %%xmm1\n\t"
				"mulps	%3, %%xmm1\n\t"
				"subps	%%xmm1, %%xmm0\n\t"
				"movaps	%%xmm3, %%xmm1\n\t"
				"mulps	%4, %%xmm1\n\t"
				"subps	%%xmm1, %%xmm0\n\t"
				"movaps	%%xmm3, %%xmm1\n\t"
				"mulps	%5, %%xmm1\n\t"
				"subps	%%xmm1, %%xmm0"
				: "=m" (vexpf_cody_waite_ln_2[0]),
				  "=m" (vexpf_cody_waite_ln_2[4]),
				  "=m" (vexpf_cody_waite_ln_2[8])
				: "m" (vexpf_cody_waite_ln_2[0]),
				  "m" (vexpf_cody_waite_ln_2[4]),
				  "m" (vexpf_cody_waite_ln_2[8])
				: "%xmm0", "%xmm1", "%xmm3");
			__asm__ __volatile__ (
				"movaps	%%xmm0, %%xmm1\n\t"
				"mulps	%7, %%xmm1\n\t"
				"addps	%8, %%xmm1\n\t"
				"mulps	%%xmm0, %%xmm1\n\t"
				"addps	%9, %%xmm1\n\t"
				"mulps	%%xmm0, %%xmm1\n\t"
				"addps	%10, %%xmm1\n\t"
				"mulps	%%xmm0, %%xmm1\n\t"
				"addps	%11, %%xmm1\n\t"
				"mulps	%%xmm0, %%xmm1\n\t"
				"addps	%12, %%xmm1\n\t"
				"mulps	%%xmm0, %%xmm1\n\t"
				"mulps	%%xmm0, %%xmm1\n\t"
				"addps	%%xmm1, %%xmm0\n\t"
				"addps	%13, %%xmm0"
				: "=m" (vexpf_minimax_coefficient[0]),
				  "=m" (vexpf_minimax_coefficient[4]),
				  "=m" (vexpf_minimax_coefficient[8]),
				  "=m" (vexpf_minimax_coefficient[12]),
				  "=m" (vexpf_minimax_coefficient[16]),
				  "=m" (vexpf_minimax_coefficient[20]),
				  "=m" (vexpf_minimax_coefficient[24])
				: "m" (*vexpf_minimax_coefficient),
				  "m" (vexpf_minimax_coefficient[4]),
				  "m" (vexpf_minimax_coefficient[8]),
				  "m" (vexpf_minimax_coefficient[12]),
				  "m" (vexpf_minimax_coefficient[16]),
				  "m" (vexpf_minimax_coefficient[20]),
				  "m" (vexpf_minimax_coefficient[24])
				: "%xmm0", "%xmm1");
#ifdef HAVE_SSE2
			__asm__ __volatile__ (
				"movdqa	%3, %%xmm1\n\t"
				"pcmpgtd	%%xmm4, %%xmm1\n\t"
				"movmskps	%%xmm1, %%eax\n\t"
				"testb	$0xf, %%al\n\t"
				"jz	0f\n\t"
				"movdqa	%4, %%xmm3\n\t"
				"pand	%%xmm1, %%xmm3\n\t"
				"paddd	%%xmm4, %%xmm3\n\t"
				"pslld	$23, %%xmm3\n\t"
				"paddd	%%xmm3, %%xmm0\n\t"
				"movaps	%%xmm0, %%xmm3\n\t"
				"mulps	%5, %%xmm3\n\t"
				"andps	%%xmm1, %%xmm3\n\t"
				"andnps	%%xmm0, %%xmm1\n\t"
				"orps	%%xmm3, %%xmm1\n\t"
				"movaps	%%xmm1, %%xmm0\n\t"
				"jmp	1f\n"
				"0:\n\t"
				"pslld	$23, %%xmm4\n\t"
				"paddd	%%xmm4, %%xmm0\n"
				"1:\n\t"
				: "=m" (*vexpf_denormal_exponent),
				  "=m" (*vexpf_denormal_exponent_shift),
				  "=m" (*vexpf_denormal_scale)
				: "m" (*vexpf_denormal_exponent),
				  "m" (*vexpf_denormal_exponent_shift),
				  "m" (*vexpf_denormal_scale)
				: "%eax", "%xmm0", "%xmm1", "%xmm3", "%xmm4");
#else // HAVE_SSE2
			__asm__ __volatile__ (
				"movaps	%%xmm0, %0\n\t"
				"pextrw	$0x0, %%xmm4, %%eax\n\t"
				"cmpw	$-126, %%ax\n\t"
				"jg	0f\n\t"
				"addw	$24, %%ax\n\t"
				"shll	$23, %%eax\n\t"
				"addl	%%eax, %0\n\t"
				"movss	%2, %%xmm0\n\t"
				"mulss	%3, %%xmm0\n\t"
				"movss	%%xmm0, %0\n"
				"0:\n\t"
				"shll	$23, %%eax\n\t"
				"addl	%%eax, %0"
				: "=m" (vexpf_buffer[0]),
				  "=m" (*vexpf_denormal_scale)
				: "m" (vexpf_buffer[0]),
				  "m" (*vexpf_denormal_scale)
				: "%eax", "%xmm0", "%xmm4");
			__asm__ __volatile__ (
				"pextrw	$0x2, %%xmm4, %%eax\n\t"
				"cmpw	$-126, %%ax\n\t"
				"jg	0f\n\t"
				"addw	$24, %%ax\n\t"
				"shll	$23, %%eax\n\t"
				"addl	%%eax, %0\n\t"
				"movss	%2, %%xmm0\n\t"
				"mulss	%3, %%xmm0\n\t"
				"movss	%%xmm0, %0\n"
				"0:\n\t"
				"shll	$23, %%eax\n\t"
				"addl	%%eax, %0"
				: "=m" (vexpf_buffer[1]),
				  "=m" (*vexpf_denormal_scale)
				: "m" (vexpf_buffer[1]),
				  "m" (*vexpf_denormal_scale)
				: "%eax", "%xmm0", "%xmm4");
			__asm__ __volatile__ (
				"pextrw	$0x4, %%xmm4, %%eax\n\t"
				"cmpw	$-126, %%ax\n\t"
				"jg	0f\n\t"
				"addw	$24, %%ax\n\t"
				"shll	$23, %%eax\n\t"
				"addl	%%eax, %0\n\t"
				"movss	%2, %%xmm0\n\t"
				"mulss	%3, %%xmm0\n\t"
				"movss	%%xmm0, %0\n"
				"0:\n\t"
				"shll	$23, %%eax\n\t"
				"addl	%%eax, %0"
				: "=m" (vexpf_buffer[2]),
				  "=m" (*vexpf_denormal_scale)
				: "m" (vexpf_buffer[2]),
				  "m" (*vexpf_denormal_scale)
				: "%eax", "%xmm0", "%xmm4");
			__asm__ __volatile__ (
				"pextrw	$0x6, %%xmm4, %%eax\n\t"
				"cmpw	$-126, %%ax\n\t"
				"jg	0f\n\t"
				"addw	$24, %%ax\n\t"
				"shll	$23, %%eax\n\t"
				"addl	%%eax, %0\n\t"
				"movss	%2, %%xmm0\n\t"
				"mulss	%3, %%xmm0\n\t"
				"movss	%%xmm0, %0\n"
				"0:\n\t"
				"shll	$23, %%eax\n\t"
				"addl	%%eax, %0"
				: "=m" (vexpf_buffer[3]),
				  "=m" (*vexpf_denormal_scale)
				: "m" (vexpf_buffer[3]),
				  "m" (*vexpf_denormal_scale)
				: "%eax", "%xmm0", "%xmm4");
#endif // HAVE_SSE2
			__asm__ __volatile__ (
				"movaps	%%xmm2, %%xmm3\n\t"
				"cmpltps	%5, %%xmm3\n\t"
#ifdef HAVE_SSE2
				"andnps	%%xmm0, %%xmm3\n\t"
#else // HAVE_SSE2
				"andnps	%7, %%xmm3\n\t"
#endif // HAVE_SSE2
				"movaps	%6, %%xmm1\n\t"
				"cmpleps	%%xmm2, %%xmm1\n\t"
				"movmskps	%%xmm1, %%eax\n\t"
				"testb	$0xf, %%al\n\t"
				"jz	0f\n\t"
				"movaps	%%xmm1, %%xmm4\n\t"
				"movaps	%8, %%xmm0\n\t"
				"addps	%%xmm0, %%xmm0\n\t"
				"andps	%%xmm0, %%xmm4\n\t"
				"andnps	%%xmm3, %%xmm1\n\t"
				"orps	%%xmm1, %%xmm4\n\t"
				"movaps	%%xmm4, %%xmm3\n"
				"0:\n\t"
				"movaps	%%xmm2, %%xmm0\n\t"
				"cmpeqps	%%xmm0, %%xmm0\n\t"
				"andps	%%xmm0, %%xmm3\n\t"
				"orps	%9, %%xmm2\n\t"
				"andnps	%%xmm2, %%xmm0\n\t"
				"orps	%%xmm3, %%xmm0"
				: "=m" (*vexpf_x_underflow),
				  "=m" (*vexpf_x_overflow),
				  "=m" (*vexpf_buffer), "=m" (*vexpf_float_max),
				  "=m" (*vexpf_qnan)
				: "m" (*vexpf_x_underflow),
				  "m" (*vexpf_x_overflow),
				  "m" (*vexpf_buffer), "m" (*vexpf_float_max),
				  "m" (*vexpf_qnan)
				: "%eax", "%xmm0", "%xmm1", "%xmm2", "%xmm3",
				  "%xmm4");
			__asm__ __volatile__ (
				"movss	%6, %%xmm1\n\t"
				"movss	%7, %%xmm2\n\t"
				"unpcklps	%%xmm2, %%xmm1\n\t"
				"movss	%8, %%xmm2\n\t"
				"movss	%9, %%xmm3\n\t"
				"unpcklps	%%xmm3, %%xmm2\n\t"
				"movlhps	%%xmm2, %%xmm1\n\t"
				"mulps	%%xmm1, %%xmm0\n\t"
				"movss	%%xmm7, %%xmm1\n\t"
				"addss	%%xmm1, %%xmm1\n\t"
				"shufps	$0x0, %%xmm1, %%xmm1\n\t"
				"movaps	%%xmm1, %%xmm3\n\t"
				"mulps	%%xmm3, %%xmm3\n\t"
				"mulps	%%xmm0, %%xmm1\n\t"
				"movaps	%%xmm1, %%xmm2\n\t"
				"mulps	%%xmm0, %%xmm3\n\t"
#ifdef __x86_64__
				"addps	%%xmm0, %%xmm10\n\t"
				"mulps	%%xmm8, %%xmm1\n\t"
				"addps	%%xmm1, %%xmm11\n\t"
				"mulps	%%xmm9, %%xmm2\n\t"
				"addps	%%xmm2, %%xmm12\n\t"
				"mulps	%%xmm8, %%xmm3\n\t"
				"mulps	%%xmm9, %%xmm3\n\t"
				"addps	%%xmm3, %%xmm13\n\t"
#else // __x86_64__
				"mulps	%10, %%xmm1\n\t"
				"mulps	%11, %%xmm2\n\t"
				"mulps	%10, %%xmm3\n\t"
				"mulps	%11, %%xmm3"
#endif // __x86_64__
				: "=m" (track_momentum_0[1]),
				  "=m" (track_momentum_1[1]),
				  "=m" (track_momentum_2[1]),
				  "=m" (track_momentum_3[1]),
				  "=m" (deviation_buffer[0]),
				  "=m" (deviation_buffer[4])
				: "m" (track_momentum_0[1]),
				  "m" (track_momentum_1[1]),
				  "m" (track_momentum_2[1]),
				  "m" (track_momentum_3[1]),
				  "m" (deviation_buffer[0]),
				  "m" (deviation_buffer[4])
				: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
				  "%xmm7"
#ifdef __x86_64__
				  , "%xmm8", "%xmm9", "%xmm10", "%xmm11", "%xmm12",
				  "%xmm13", "%xmm14", "%xmm15"
#endif // __x86_64__
				);
#ifndef __x86_64__
			__asm__ __volatile__ (
#ifdef HAVE_SSE3
				"haddps	%%xmm1, %%xmm1\n\t"
				"haddps	%%xmm1, %%xmm1\n\t"
				"haddps	%%xmm2, %%xmm2\n\t"
				"haddps	%%xmm2, %%xmm2\n\t"
				"haddps	%%xmm3, %%xmm3\n\t"
				"haddps	%%xmm3, %%xmm3\n\t"
#else // HAVE_SSE3
				"movhlps	%%xmm1, %%xmm4\n\t"
				"addps	%%xmm1, %%xmm4\n\t"
				"movaps	%%xmm4, %%xmm1\n\t"
				"shufps	$0x55, %%xmm1, %%xmm1\n\t"
				"addss	%%xmm4, %%xmm1\n\t"
				"movhlps	%%xmm2, %%xmm4\n\t"
				"addps	%%xmm2, %%xmm4\n\t"
				"movaps	%%xmm4, %%xmm2\n\t"
				"shufps	$0x55, %%xmm2, %%xmm2\n\t"
				"addss	%%xmm4, %%xmm2\n\t"
				"movhlps	%%xmm3, %%xmm4\n\t"
				"addps	%%xmm3, %%xmm4\n\t"
				"movaps	%%xmm4, %%xmm3\n\t"
				"shufps	$0x55, %%xmm3, %%xmm3\n\t"
				"addss	%%xmm4, %%xmm3\n\t"
#endif // HAVE_SSE3
				"unpcklps	%%xmm2, %%xmm1\n\t"
				"movlhps	%%xmm1, %%xmm4\n\t"
#ifdef HAVE_SSE3
				"movaps	%%xmm0, %%xmm1\n\t"
				"haddps	%%xmm1, %%xmm1\n\t"
				"haddps	%%xmm1, %%xmm1\n\t"
#else // HAVE_SSE3
				"movhlps	%%xmm0, %%xmm1\n\t"
				"addps	%%xmm0, %%xmm1\n\t"
				"movaps	%%xmm1, %%xmm2\n\t"
				"shufps	$0x55, %%xmm2, %%xmm2\n\t"
				"addss	%%xmm2, %%xmm1\n\t"
#endif // HAVE_SSE3
				"movss	%%xmm1, %%xmm4\n\t"
				"addps	%%xmm4, %%xmm5"
				: : : "%xmm0", "%xmm1", "%xmm2", "%xmm4", "%xmm5");
#endif // __x86_64__

			static float one[4] __attribute__ ((aligned(16))) = {
				1.0F, 1.0F, 1.0F, 1.0F
			};

#ifdef __x86_64__
			__asm__ __volatile__ (
				"movss	%%xmm7, %%xmm1\n\t"
				"addss	%%xmm1, %%xmm1\n\t"
				"shufps	$0x0, %%xmm1, %%xmm1\n\t"
				"mulps	%%xmm8, %%xmm8\n\t"
				"mulps	%%xmm1, %%xmm8\n\t"
				"addps	%1, %%xmm8\n\t"
				"mulps	%%xmm1, %%xmm8\n\t"
				"mulps	%%xmm0, %%xmm8\n\t"
				"addps	%%xmm8, %%xmm14\n\t"
				"mulps	%%xmm9, %%xmm9\n\t"
				"mulps	%%xmm1, %%xmm9\n\t"
				"addps	%1, %%xmm9\n\t"
				"mulps	%%xmm1, %%xmm9\n\t"
				"mulps	%%xmm0, %%xmm9\n\t"
				"addps	%%xmm9, %%xmm15\n\t"
				: "=m" (*one)
				: "m" (*one)
				: "%xmm0", "%xmm8", "%xmm9", "%xmm14", "%xmm15");
#else // __x86_64__
			__asm__ __volatile__ (
				"movss	%%xmm7, %%xmm1\n\t"
				"addss	%%xmm1, %%xmm1\n\t"
				"shufps	$0x0, %%xmm1, %%xmm1\n\t"
#ifdef __x86_64__
				"mulps	%%xmm8, %%xmm8\n\t"
				"mulps	%%xmm1, %%xmm8\n\t"
				"addps	%3, %%xmm8\n\t"
				"mulps	%%xmm1, %%xmm8\n\t"
				"mulps	%%xmm0, %%xmm8\n\t"
				"mulps	%%xmm9, %%xmm9\n\t"
				"mulps	%%xmm1, %%xmm9\n\t"
				"addps	%3, %%xmm9\n\t"
				"mulps	%%xmm1, %%xmm9\n\t"
				"mulps	%%xmm0, %%xmm9\n\t"
#ifdef HAVE_SSE3
				"haddps	%%xmm8, %%xmm8\n\t"
				"haddps	%%xmm8, %%xmm8\n\t"
				"haddps	%%xmm9, %%xmm9\n\t"
				"haddps	%%xmm9, %%xmm9\n\t"
#else // HAVE_SSE3
				"movhlps	%%xmm8, %%xmm1\n\t"
				"addps	%%xmm8, %%xmm1\n\t"
				"movaps	%%xmm1, %%xmm8\n\t"
				"shufps	$0x55, %%xmm8, %%xmm8\n\t"
				"addss	%%xmm1, %%xmm8\n\t"
				"movhlps	%%xmm9, %%xmm1\n\t"
				"addps	%%xmm9, %%xmm1\n\t"
				"movaps	%%xmm1, %%xmm9\n\t"
				"shufps	$0x55, %%xmm9, %%xmm9\n\t"
				"addss	%%xmm1, %%xmm9\n\t"
#endif // HAVE_SSE3
				"unpcklps	%%xmm9, %%xmm8\n\t"
				"movlhps	%%xmm3, %%xmm8\n\t"
				"addps	%%xmm8, %%xmm6"
#else // __x86_64__
				"movaps	%4, %%xmm2\n\t"
				"mulps	%%xmm2, %%xmm2\n\t"
				"mulps	%%xmm1, %%xmm2\n\t"
				"addps	%3, %%xmm2\n\t"
				"mulps	%%xmm1, %%xmm2\n\t"
				"mulps	%%xmm0, %%xmm2\n\t"
				"movaps	%5, %%xmm4\n\t"
				"mulps	%%xmm4, %%xmm4\n\t"
				"mulps	%%xmm1, %%xmm4\n\t"
				"addps	%3, %%xmm4\n\t"
				"mulps	%%xmm1, %%xmm4\n\t"
				"mulps	%%xmm0, %%xmm4\n\t"
#ifdef HAVE_SSE3
				"haddps	%%xmm2, %%xmm2\n\t"
				"haddps	%%xmm2, %%xmm2\n\t"
				"haddps	%%xmm4, %%xmm4\n\t"
				"haddps	%%xmm4, %%xmm4\n\t"
#else // HAVE_SSE3
				"movhlps	%%xmm2, %%xmm1\n\t"
				"addps	%%xmm2, %%xmm1\n\t"
				"movaps	%%xmm1, %%xmm2\n\t"
				"shufps	$0x55, %%xmm2, %%xmm2\n\t"
				"addss	%%xmm1, %%xmm2\n\t"
				"movhlps	%%xmm4, %%xmm1\n\t"
				"addps	%%xmm4, %%xmm1\n\t"
				"movaps	%%xmm1, %%xmm4\n\t"
				"shufps	$0x55, %%xmm4, %%xmm4\n\t"
				"addss	%%xmm1, %%xmm4\n\t"
#endif // HAVE_SSE3
				"unpcklps	%%xmm4, %%xmm2\n\t"
				"movlhps	%%xmm3, %%xmm2\n\t"
				"addps	%%xmm2, %%xmm6"
#endif // __x86_64__
				: "=m" (*one), "=m" (deviation_buffer[0]),
				  "=m" (deviation_buffer[4])
				: "m" (*one), "m" (deviation_buffer[0]),
				  "m" (deviation_buffer[4])
				: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
				  "%xmm6", "%xmm7"
#ifdef __x86_64__
				  , "%xmm8", "%xmm9"
#endif // __x86_64__
				);
#endif // __x86_64__
		}

#ifdef __x86_64__
		__asm__ __volatile__ (
#ifdef HAVE_SSE3
			// SSE3 (in place) horizontal add
			"haddps	%%xmm10, %%xmm10\n\t"
			"haddps	%%xmm10, %%xmm10\n\t"
			// %xmm10 -> %xmm11
			"haddps	%%xmm11, %%xmm11\n\t"
			"haddps	%%xmm11, %%xmm11\n\t"
			// %xmm11 -> %xmm12
			"haddps	%%xmm12, %%xmm12\n\t"
			"haddps	%%xmm12, %%xmm12\n\t"
			// %xmm12 -> %xmm13
			"haddps	%%xmm13, %%xmm13\n\t"
			"haddps	%%xmm13, %%xmm13\n\t"
			// %xmm13 -> %xmm14
			"haddps	%%xmm14, %%xmm14\n\t"
			"haddps	%%xmm14, %%xmm14\n\t"
			// %xmm14 -> %xmm15
			"haddps	%%xmm15, %%xmm15\n\t"
			"haddps	%%xmm15, %%xmm15\n\t"
#else // HAVE_SSE3
			// SSE (in place, simulated) horizontal add
			"movhlps	%%xmm10, %%xmm0\n\t"
			"addps	%%xmm10, %%xmm0\n\t"
			"movaps	%%xmm0, %%xmm10\n\t"
			"shufps	$0x55, %%xmm10, %%xmm10\n\t"
			"addss	%%xmm0, %%xmm10\n\t"
			// %xmm10 -> %xmm11
			"movhlps	%%xmm11, %%xmm0\n\t"
			"addps	%%xmm11, %%xmm0\n\t"
			"movaps	%%xmm0, %%xmm11\n\t"
			"shufps	$0x55, %%xmm11, %%xmm11\n\t"
			"addss	%%xmm0, %%xmm11\n\t"
			// %xmm11 -> %xmm12
			"movhlps	%%xmm12, %%xmm0\n\t"
			"addps	%%xmm12, %%xmm0\n\t"
			"movaps	%%xmm0, %%xmm12\n\t"
			"shufps	$0x55, %%xmm12, %%xmm12\n\t"
			"addss	%%xmm0, %%xmm12\n\t"
			// %xmm12 -> %xmm13
			"movhlps	%%xmm13, %%xmm0\n\t"
			"addps	%%xmm13, %%xmm0\n\t"
			"movaps	%%xmm0, %%xmm13\n\t"
			"shufps	$0x55, %%xmm13, %%xmm13\n\t"
			"addss	%%xmm0, %%xmm13\n\t"
			// %xmm13 -> %xmm14
			"movhlps	%%xmm14, %%xmm0\n\t"
			"addps	%%xmm14, %%xmm0\n\t"
			"movaps	%%xmm0, %%xmm14\n\t"
			"shufps	$0x55, %%xmm14, %%xmm14\n\t"
			"addss	%%xmm0, %%xmm14\n\t"
			// %xmm14 -> %xmm15
			"movhlps	%%xmm15, %%xmm0\n\t"
			"addps	%%xmm15, %%xmm0\n\t"
			"movaps	%%xmm0, %%xmm15\n\t"
			"shufps	$0x55, %%xmm15, %%xmm15\n\t"
			"addss	%%xmm0, %%xmm15\n\t"
#endif // HAVE_SSE3
			"movss	%%xmm10, %0\n\t"
			"unpcklps	%%xmm12, %%xmm11\n\t"
			"movlps	%%xmm11, %1\n\t"
			"unpcklps	%%xmm15, %%xmm14\n\t"
			"movlhps	%%xmm13, %%xmm14\n\t"
			"movaps	%%xmm14, %2"
			: "=m" (perp), "=m" (*gradient), "=m" (*hessian)
			: "m" (perp), "m" (*gradient), "m" (*hessian)
			: "%xmm10", "%xmm11", "%xmm12", "%xmm13", "%xmm14",
			  "%xmm15");
#else // __x86_64__
		__asm__ __volatile__ (
			"movss	%%xmm5, %0\n\t"
			"movhps	%%xmm5, %1\n\t"
			"movaps	%%xmm6, %2"
			: "=m" (perp), "=m" (*gradient), "=m" (*hessian)
			: "m" (perp), "m" (*gradient), "m" (*hessian)
			: "%xmm5", "%xmm6");
#endif // __x86_64__

#ifndef NVERIFY
		ID(float perp_verify = 0);
		ID(float gradient_verify[2]);
		IS(gradient_verify[0] = 0);
		IS(gradient_verify[1] = 0);
		ID(float hessian_verify[3]);
		IS(hessian_verify[0] = 0);
		IS(hessian_verify[1] = 0);
		IS(hessian_verify[2] = 0);
		IS(iterator = track_begin);
		for(int i = track_end - track_begin; i > 0; i--) {
			ID(const float deviation_pseudorapidity =
			   jet.momentum().pseudorapidity() -
			   iterator->momentum().pseudorapidity());
			ID(const float deviation_azimuth =
			   angular_range_reduce(jet.momentum().azimuth() -
									iterator->momentum().azimuth()));
			ID(const float radial_deviation =
			   deviation_pseudorapidity *
			   deviation_pseudorapidity +
			   deviation_azimuth * deviation_azimuth);
			ID(const float exponential_factor =
			   iterator->momentum().perp() *
			   exp(scale * radial_deviation));
			IS(perp_verify += exponential_factor);
			IS(gradient_verify[0] +=
			   2.0F * scale * deviation_pseudorapidity *
			   exponential_factor);
			IS(gradient_verify[1] +=
			   2.0F * scale * deviation_azimuth *
			   exponential_factor);
			IS(hessian_verify[0] += 2.0F * scale *
			   (1.0F + 2.0F * scale * deviation_pseudorapidity *
				deviation_pseudorapidity) *
			   exponential_factor);
			IS(hessian_verify[2] += 4.0F * scale * scale *
			   deviation_azimuth * deviation_pseudorapidity *
			   exponential_factor);
			IS(hessian_verify[1] += 2.0F * scale *
			   (1.0F + 2.0F * scale * deviation_azimuth *
				deviation_azimuth) * exponential_factor);
			IS(iterator++);
		}
		I(FEQ(perp, perp_verify));
		I(A(int i = 0, i < 2, i++,
			FEQ(gradient[i], gradient_verify[i])));
		I(A(int i = 0, i < 3, i++,
			FEQ(hessian[i], hessian_verify[i])));
#endif // NVERIFY

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

#ifdef __INTEL_COMPILER
#pragma warning(disable: 869)
#endif // __INTEL_COMPILER

#endif // HAVE_SSE
}
