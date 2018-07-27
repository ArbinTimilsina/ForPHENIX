#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H

#define _XOPEN_SOURCE 600
#include <cstdlib>
#include <jetrec/rec.h>
#include <jetrec/util.h>

namespace jet {

#ifdef HAVE_SSE

	// SSE/SSE2/SSE3 implementation of filtering operations for
	// IA-32/AMD64/Intel 64 architectures

	void reconstruction_filtering_t::
	maximum_map_sse(void) const
	{
		const int _stride = _npixel_azimuth;
		const int vector_length = 4;

		I(A(int i = _stride, i < _npixel_pseudorapidity * _stride,
			i++, F(_distribution[i])));

		// Convolution based SSE maximum finder

		const int npixel_pseudorapidity_mod_3 =
			_npixel_pseudorapidity % 3;
		const int unrolled_count = (_npixel_pseudorapidity - 2) / 3;

		for(int row_index = 0; row_index < _stride;
			row_index += vector_length) {
			float *pmaximum_2 = _maximum_2 + row_index;
			float *pmaximum_3 = _maximum_3 + row_index;

			// FIXME: Pseudorapidity boundary values must be properly
			// initialized (they are not used for now).

			__asm__ __volatile__ (
				"xorps	%%xmm1, %%xmm1\n\t"
				"movaps	%%xmm1, %0\n\t"
				"movaps	%%xmm1, %1"
				: "=m" (*pmaximum_2), "=m" (*pmaximum_3)
				: "m" (*pmaximum_2), "m" (*pmaximum_3)
				: "%xmm1");
			pmaximum_2 += _stride;
			pmaximum_3 += _stride;

			float *pdistribution = _distribution + row_index;

			// Initialize the ring buffer minus the
			// _distribution[j:j+vl-1:1,2] which is going to be loaded
			// inside the loop.
			__asm__ __volatile__ (
				// The common index (j : j + vl - 1 : 1) is suppressed
				// below.
				// Load _distribution[0]
				"movaps	%1, %%xmm1"
				: "=m" (*pdistribution)
				: "m" (*pdistribution)
				: "%xmm1", "%xmm2");
			pdistribution += _stride;
			__asm__ __volatile__ (
				// Load _distribution[1]
				"movaps	%1, %%xmm2"
				: "=m" (*pdistribution)
				: "m" (*pdistribution)
				: "%xmm1", "%xmm2");
			pdistribution += _stride;

			// FIXME: Are there speed to be gain by fitting as much
			// stripes as there are SSE registers (2 for IA-32 and 5
			// for AMD64/Intel 64) and improving memory locality?

			for(int i = unrolled_count; i > 0; i--) {
				// Use %xmm0, %xmm1, and %xmm2 as a register ring
				// buffer and unroll the loop to the ring buffer size,
				// so that the register permutations are loop
				// invariant.
				__asm__ __volatile__ (
					// The common index (j : j + vl - 1 : 1) is
					// suppressed below.
					// Load _distribution[i]
					"movaps	%3, %%xmm0\n\t"
					// %xmm0 = _distribution[i]
					// %xmm1 = _distribution[i - 2]
					// %xmm2 = _distribution[i - 1]
					// %xmm1 = max(%xmm0, %xmm1)
					"maxps	%%xmm0, %%xmm1\n\t"
					// Store _maximum_2[i - 1]
					"movaps	%%xmm1, %1\n\t"
					// %xmm1 = max(%xmm2, %xmm1)
					"maxps	%%xmm2, %%xmm1\n\t"
					// Store _maximum_3[i - 1]
					"movaps	%%xmm1, %2"
					: "=m" (*pdistribution), "=m" (*pmaximum_2),
					  "=m" (*pmaximum_3)
					: "m" (*pdistribution), "m" (*pmaximum_2),
					  "m" (*pmaximum_3)
					: "%xmm0", "%xmm1", "%xmm2");
				pdistribution += _stride;
				pmaximum_2 += _stride;
				pmaximum_3 += _stride;
				__asm__ __volatile__ (
					// Cyclically permute %xmm0 -> %xmm1 -> %xmm2 ->
					// %xmm0
					"movaps	%3, %%xmm1\n\t"
					"maxps	%%xmm1, %%xmm2\n\t"
					"movaps	%%xmm2, %1\n\t"
					"maxps	%%xmm0, %%xmm2\n\t"
					"movaps	%%xmm2, %2"
					: "=m" (*pdistribution), "=m" (*pmaximum_2),
					  "=m" (*pmaximum_3)
					: "m" (*pdistribution), "m" (*pmaximum_2),
					  "m" (*pmaximum_3)
					: "%xmm0", "%xmm1", "%xmm2");
				pdistribution += _stride;
				pmaximum_2 += _stride;
				pmaximum_3 += _stride;
				__asm__ __volatile__ (
					// Cyclically permute %xmm0 -> %xmm1 -> %xmm2 ->
					// %xmm0
					"movaps	%3, %%xmm2\n\t"
					"maxps	%%xmm2, %%xmm0\n\t"
					"movaps	%%xmm0, %1\n\t"
					"maxps	%%xmm1, %%xmm0\n\t"
					"movaps	%%xmm0, %2"
					: "=m" (*pdistribution), "=m" (*pmaximum_2),
					  "=m" (*pmaximum_3)
					: "m" (*pdistribution), "m" (*pmaximum_2),
					  "m" (*pmaximum_3)
					: "%xmm0", "%xmm1", "%xmm2");
				pdistribution += _stride;
				pmaximum_2 += _stride;
				pmaximum_3 += _stride;
			}
			// Calculate the remaining segment.
			switch(npixel_pseudorapidity_mod_3) {
			case 0:
				__asm__ __volatile__ (
					// Part 1 from the loop above, fusing one load
					// with maxps
					"maxps	%3, %%xmm1\n\t"
					"movaps	%%xmm1, %1\n\t"
					"maxps	%%xmm2, %%xmm1\n\t"
					"movaps	%%xmm1, %2"
					: "=m" (*pdistribution), "=m" (*pmaximum_2),
					  "=m" (*pmaximum_3)
					: "m" (*pdistribution), "m" (*pmaximum_2),
					  "m" (*pmaximum_3)
					: "%xmm0", "%xmm1", "%xmm2");
				pmaximum_2 += _stride;
				pmaximum_3 += _stride;
				break;
			case 1:
				__asm__ __volatile__ (
					// Part 1 from the loop above
					"movaps	%3, %%xmm0\n\t"
					"maxps	%%xmm0, %%xmm1\n\t"
					"movaps	%%xmm1, %1\n\t"
					"maxps	%%xmm2, %%xmm1\n\t"
					"movaps	%%xmm1, %2"
					: "=m" (*pdistribution), "=m" (*pmaximum_2),
					  "=m" (*pmaximum_3)
					: "m" (*pdistribution), "m" (*pmaximum_2),
					  "m" (*pmaximum_3)
					: "%xmm0", "%xmm1", "%xmm2");
				pdistribution += _stride;
				pmaximum_2 += _stride;
				pmaximum_3 += _stride;
				__asm__ __volatile__ (
					// Part 2 from the loop above, fusing one load
					// with maxps
					"maxps	%3, %%xmm2\n\t"
					"movaps	%%xmm2, %1\n\t"
					"maxps	%%xmm0, %%xmm2\n\t"
					"movaps	%%xmm2, %2"
					: "=m" (*pdistribution), "=m" (*pmaximum_2),
					  "=m" (*pmaximum_3)
					: "m" (*pdistribution), "m" (*pmaximum_2),
					  "m" (*pmaximum_3)
					: "%xmm0", "%xmm1", "%xmm2");
				pmaximum_2 += _stride;
				pmaximum_3 += _stride;
				break;
			}
			__asm__ __volatile__ (
				"xorps	%%xmm1, %%xmm1\n\t"
				"movaps	%%xmm1, %0\n\t"
				"movaps	%%xmm1, %1"
				: "=m" (*pmaximum_2), "=m" (*pmaximum_3)
				: "m" (*pmaximum_2), "m" (*pmaximum_3)
				: "%xmm1");
		}

#ifndef NVERIFY
		// Postcondition for the inline assembly code above, quite
		// expensive to recalculate.
		I(A(int i = _stride,
			i < (_npixel_pseudorapidity - 1) * _stride, i++,
			_maximum_2[i] == std::max(_distribution[i - _stride],
									  _distribution[i + _stride]) &&
			_maximum_3[i] == std::max(_distribution[i],
									  _maximum_2[i])));
#endif // NVERIFY

		int index = 0;

		for(int row_index = 0; row_index < _npixel_azimuth;
			row_index++) {
			_maximum[index] = false;
			index++;
		}
		for(int column_index = 1;
			column_index < _npixel_pseudorapidity - 1;
			column_index++) {
			I(index + _npixel_azimuth - 1 >= 0 &&
			  index + _npixel_azimuth - 1 < _npixel);
			I(index >= 0 && index < _npixel);
			I(index + 1 >= 0 && index + 1 < _npixel);

			_maximum[index] =
				_distribution[index] >
				_maximum_3[index + _npixel_azimuth - 1] &&
				_distribution[index] > _maximum_2[index] &&
				_distribution[index] > _maximum_3[index + 1];
			index++;
			for(int row_index = 1; row_index < _npixel_azimuth - 1;
				row_index++) {
				I(index == _stride * column_index + row_index);
				I(index - 1 >= 0 && index - 1 < _npixel);
				I(index >= 0 && index < _npixel);
				I(index + 1 >= 0 && index + 1 < _npixel);

				_maximum[index] =
					_distribution[index] > _maximum_3[index - 1] &&
					_distribution[index] > _maximum_2[index] &&
					_distribution[index] > _maximum_3[index + 1];
				index++;
			}

			I(index - 1 >= 0 && index - 1 < _npixel);
			I(index >= 0 && index < _npixel);
			I(index - _npixel_azimuth + 1 >= 0 &&
			  index - _npixel_azimuth + 1 < _npixel);

			_maximum[index] =
				_distribution[index] > _maximum_3[index - 1] &&
				_distribution[index] > _maximum_2[index] &&
				_distribution[index] >
				_maximum_3[index - _npixel_azimuth + 1];
			index++;
		}
		for(int row_index = 0; row_index < _npixel_azimuth;
			row_index++) {
			_maximum[index] = false;
			index++;
		}
	}

#endif // HAVE_SSE

#ifdef HAVE_SSE2

#if 0
	// INCOMPLETE, DO NOT USE!
	void reconstruction_filtering_t::
	stationary_map_sse2(void) const
	{
		const int _stride = _npixel_azimuth;
		const int vector_length = 4;

		I(A(int i = _stride, i < _npixel_pseudorapidity * _stride,
			i++, F(_distribution[i])));

		for(int row_index = 0; row_index < _stride;
			row_index += vector_length) {
			bool *pmaximum = _maximum + row_index + _stride;
			float *pdistribution = _distribution + row_index;

			__asm__ __volatile__ (
				"movaps	%1, %%xmm0"
				: "=m" (*pdistribution)
				: "m" (*pdistribution)
				: "%xmm0");
			pdistribution += _stride;

			// FIXME: Are there speed to be gain by fitting as much
			// stripes as there are SSE registers (2 for IA-32 and 5
			// for AMD64/Intel 64) and improving memory locality?

			for(int i = (_npixel_pseudorapidity - 1) >> 2; i > 0;
				i--) {
				__asm__ __volatile__ (
					"movaps	%2, %%xmm1\n\t"
					"cmpeqps	%%xmm1, %%xmm0\n\t"
					// Unpack and store %xmm0
					"packusdw	%%xmm0, %%xmm0\n\t"
					"packuswb	%%xmm0, %%xmm0\n\t"
					"movq	%3, %%xmm2\n\t"
					"por	%%xmm2, %%xmm0\n\t"
					"movq	%%xmm0, %1"
					: "=m" (*pdistribution), "=m" (*pmaximum)
					: "m" (*pdistribution), "m" (*pmaximum)
					: "%xmm0", "%xmm1", "%xmm2");
				pdistribution += _stride;
				pmaximum += _stride;
				__asm__ __volatile__ (
					// %xmm0 <-> %xmm1
					"movaps	%2, %%xmm0\n\t"
					"cmpeqps	%%xmm0, %%xmm1\n\t"
					"packusdw	%%xmm1, %%xmm1\n\t"
					"packuswb	%%xmm1, %%xmm1\n\t"
					"movq	%3, %%xmm2\n\t"
					"por	%%xmm2, %%xmm1\n\t"
					"movq	%%xmm1, %1"
					: "=m" (*pdistribution), "=m" (*pmaximum)
					: "m" (*pdistribution), "m" (*pmaximum)
					: "%xmm0", "%xmm1", "%xmm2");
				pdistribution += _stride;
				pmaximum += _stride;
			}
			// Calculate the remaining segment.
			if(_npixel_pseudorapidity & 1 == 0) {
				__asm__ __volatile__ (
					"movaps	%2, %%xmm1\n\t"
					"cmpeqps	%%xmm1, %%xmm0\n\t"
					// Unpack and store %xmm0
					"packusdw	%%xmm0, %%xmm0\n\t"
					"packuswb	%%xmm0, %%xmm0\n\t"
					"movq	%3, %%xmm2\n\t"
					"por	%%xmm2, %%xmm0\n\t"
					"movq	%%xmm0, %1"
					: "=m" (*pdistribution), "=m" (*pmaximum)
					: "m" (*pdistribution), "m" (*pmaximum)
					: "%xmm0", "%xmm1", "%xmm2");
			}
		}
	}
#endif

#endif // HAVE_SSE2

}
