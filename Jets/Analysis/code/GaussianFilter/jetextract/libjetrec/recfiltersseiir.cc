#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H

#define _XOPEN_SOURCE 600
#include <cstdlib>
#include <jetrec/rec.h>
#include <jetrec/util.h>

namespace jet {

#ifdef HAVE_SSE

	void reconstruction_filtering_iir_t::
	apply_forward_aos_2(float output[], float input[], const int n,
						float *coefficient_1, float *coefficient_2)
	{
		// "Fused" AOS forward filtering

		I(A(int i = 0, i < n, i++, F(input[i])));
		I(A(int i = 0, i < 5, i++, F(coefficient_1[i])));
		I(A(int i = 0, i < 5, i++, F(coefficient_2[i])));

		if(std::fpclassify(coefficient_1[0]) == FP_ZERO &&
		   std::fpclassify(coefficient_2[0]) == FP_ZERO) {
			// While declaring a shuffled array like the following
			// (then movaps _vector_coefficient[0], %xmm4; movaps
			// _vector_coefficient[4], %xmm5) might work for compilers
			// like GCC, Intel C++ Compiler for example will optimize
			// it away, and pass the original pointer to movaps (which
			// most likely will end in a segfault):
			//
			// float vector_coefficient[8]
			// 	__attribute__ ((aligned(16))) = {
			// 	coefficient_1[2], coefficient_2[2],
			// 	coefficient_1[4], coefficient_2[4],
			// 	coefficient_1[1], coefficient_2[1],
			// 	coefficient_1[3], coefficient_2[3]
			// };
			__asm__ __volatile__ (
				// Load the coefficients
				"movss	%8, %%xmm4\n\t"
				// Now %xmm4 = [ * | * | * | coefficient_1[2] ]
				"movss	%9, %%xmm3\n\t"
				// Now %xmm3 = [ * | * | * | coefficient_1[4] ]
				"unpcklps	%%xmm3, %%xmm4\n\t"
				// Now %xmm4 = [ * | * | coefficient_1[4] |
				// coefficient_1[2] ]
				// %xmm3 is now free
				// Same as above with %xmm3 -> %xmm1, %xmm4 -> %xmm3,
				// and for coefficient_2
				"movss	%10, %%xmm3\n\t"
				"movss	%11, %%xmm1\n\t"
				"unpcklps	%%xmm1, %%xmm3\n\t"
				// Now %xmm3 = [ * | * | coefficient_2[4] |
				// coefficient_2[2] ]
				"unpcklps	%%xmm3, %%xmm4\n\t"
				// Now %xmm4 = [ coefficient_2[3] | coefficient_1[3] |
				// coefficient_2[0] | coefficient_1[0] ]
				// Same as above, with %xmm4 <-> %xmm5 and loading
				// elements 1, 3
				"movss	%12, %%xmm5\n\t"
				"movss	%13, %%xmm3\n\t"
				"unpcklps	%%xmm3, %%xmm5\n\t"
				"movss	%14, %%xmm3\n\t"
				"movss	%15, %%xmm1\n\t"
				"unpcklps	%%xmm1, %%xmm3\n\t"
				"unpcklps	%%xmm3, %%xmm5\n\t"
				// Zero %xmm0, %xmm1
				"xorps	%%xmm0, %%xmm0\n\t"
				"xorps	%%xmm1, %%xmm1"
				: "=m" (coefficient_1[2]), "=m" (coefficient_1[4]),
				  "=m" (coefficient_2[2]), "=m" (coefficient_2[4]),
				  "=m" (coefficient_1[1]), "=m" (coefficient_1[3]),
				  "=m" (coefficient_2[1]), "=m" (coefficient_2[3])
				: "m" (coefficient_1[2]), "m" (coefficient_1[4]),
				  "m" (coefficient_2[2]), "m" (coefficient_2[4]),
				  "m" (coefficient_1[1]), "m" (coefficient_1[3]),
				  "m" (coefficient_2[1]), "m" (coefficient_2[3])
				: "%xmm0", "%xmm1", "%xmm3", "%xmm4", "%xmm5");

			float *pinput = input;

			for(int i = n >> 1; i > 0; i--) {
				__asm__ __volatile__ (
					// Start with:
					// %xmm0 = [ y2[i - 2] | y1[i - 2] | x[i - 2] |
					// x[i - 2] ]
					// %xmm1 = [ y2[i - 1] | y1[i - 1] | x[i - 1] |
					// x[i - 1] ]
					"movaps	%%xmm0, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					// Now:
					// %xmm2 = [ c2[4] * y2[i - 2] | c1[4] * y1[i - 2]
					// | c2[2] * x[i - 2] | c1[2] * x[i - 2] ]
					"movaps	%%xmm1, %%xmm3\n\t"
					"mulps	%%xmm5, %%xmm3\n\t"
					// Now:
					// %xmm3 = [ c2[3] * y2[i - 1] | c1[3] * y1[i - 1]
					// | c2[1] * x[i - 1] | c1[1] * x[i - 1] ]
					"addps	%%xmm2, %%xmm3\n\t"
					// Now:
					// %xmm3 = [ c2[3] * y2[i - 1] + c2[4] * y2[i - 2]
					// | c1[3] * y1[i - 1] + c1[4] * y1[i - 2] | c2[1]
					// * x[i - 1] + c2[2] * x[i - 2] | c1[1] * x[i -
					// 1] + c1[2] * x[i - 2] ]
					// When adapting for SSE3, note that this is not
					// an emulation of hsubps! There is also not much
					// advantage in rearranging the register layout to
					// accomodate hsubps, since the interleaved format
					// would result in another 6 cycles of latency for
					// a second shufps when loading x[i], offsetting
					// the 6 cycles of movhlps.
					"movhlps	%%xmm3, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm3\n\t"
					// Now:
					// %xmm3 = [ * | * | y2[i] | y1[i] ]
					// Load x[i]
					"movss	%1, %%xmm0\n\t"
					// %xmm0 = [ * | * | * | x[i] ]
					"shufps	$0x40, %%xmm3, %%xmm0"
					// %xmm0 = [ y2[i] | y1[i] | x[i] | x[i] ]
					: "=m" (*pinput)
					: "m" (*pinput)
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5");
				pinput++;
				__asm__ __volatile__ (
					// %xmm0 <-> %xmm1
					"movaps	%%xmm1, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"movaps	%%xmm0, %%xmm3\n\t"
					"mulps	%%xmm5, %%xmm3\n\t"
					"addps	%%xmm2, %%xmm3\n\t"
					"movhlps	%%xmm3, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm3\n\t"
					"movss	%1, %%xmm1\n\t"
					"shufps	$0x40, %%xmm3, %%xmm1"
					: "=m" (*pinput)
					: "m" (*pinput)
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5");
				pinput++;
			}
			pinput = input;

			float *poutput = output;

			for(int i = n >> 1; i > 0; i--) {
				__asm__ __volatile__ (
					// Start analogously to above
					"movaps	%%xmm0, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"movaps	%%xmm1, %%xmm3\n\t"
					"mulps	%%xmm5, %%xmm3\n\t"
					"addps	%%xmm2, %%xmm3\n\t"
					"movhlps	%%xmm3, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm3\n\t"
					// Now:
					// %xmm3 = [ * | * | y2[i] | y1[i] ]
					"movaps	%%xmm3, %%xmm2\n\t"
#ifdef HAVE_SSE3
					"haddps	%%xmm2, %%xmm2\n\t"
#else // HAVE_SSE3
					// Emulate haddps %xmm2, %xmm2
					"shufps	$0x55, %%xmm2, %%xmm2\n\t"
					// %xmm2 = [ * | * | * | y2[i] ]
					"addss	%%xmm3, %%xmm2\n\t"
#endif // HAVE_SSE3
					// %xmm2 = [ * | * | * | y1[i] + y2[i] ]
					"movss	%%xmm2, %0\n\t"
					// Load x[i]
					"movss	%3, %%xmm0\n\t"
					// %xmm0 = [ * | * | * | x[i] ]
					"shufps	$0x40, %%xmm3, %%xmm0"
					// %xmm0 = [ y2[i] | y1[i] | x[i] | x[i] ]
					// Note that the input/output are flipped
					: "=m" (*poutput), "=m" (*pinput)
					: "m" (*poutput), "m" (*pinput)
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5", "memory");
				pinput++;
				poutput++;
				__asm__ __volatile__ (
					// %xmm0 <-> %xmm1
					"movaps	%%xmm1, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"movaps	%%xmm0, %%xmm3\n\t"
					"mulps	%%xmm5, %%xmm3\n\t"
					"addps	%%xmm2, %%xmm3\n\t"
					"movhlps	%%xmm3, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm3\n\t"
					"movaps	%%xmm3, %%xmm2\n\t"
#ifdef HAVE_SSE3
					"haddps	%%xmm2, %%xmm2\n\t"
#else // HAVE_SSE3
					"shufps	$0x55, %%xmm2, %%xmm2\n\t"
					"addss	%%xmm3, %%xmm2\n\t"
#endif // HAVE_SSE3
					"movss	%%xmm2, %0\n\t"
					"movss	%3, %%xmm1\n\t"
					"shufps	$0x40, %%xmm3, %%xmm1"
					// Note that the input/output are flipped
					: "=m" (*poutput), "=m" (*pinput)
					: "m" (*poutput), "m" (*pinput)
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5", "memory");
				pinput++;
				poutput++;
			}
		}
		else {
			// assuming coefficient_1[2] == 0 && coefficient_2[2] == 0

			// Kept for reference and illustration of the shuffling
			// result, see note above.
			//
			// float vector_coefficient[8]
			// 	__attribute__ ((aligned(16))) = {
			// 	coefficient_1[0], coefficient_2[0],
			// 	coefficient_1[3], coefficient_2[3],
			// 	coefficient_1[1], coefficient_2[1],
			// 	coefficient_1[4], coefficient_2[4]
			// };
			__asm__ __volatile__ (
				// Load the coefficients
				"movss	%8, %%xmm4\n\t"
				// Now %xmm4 = [ * | * | * | coefficient_1[0] ]
				"movss	%9, %%xmm3\n\t"
				// Now %xmm3 = [ * | * | * | coefficient_1[3] ]
				"unpcklps	%%xmm3, %%xmm4\n\t"
				// Now %xmm4 = [ * | * | coefficient_1[3] |
				// coefficient_1[0] ]
				// %xmm3 is now free
				// Same as above with %xmm3 -> %xmm1, %xmm4 -> %xmm3,
				// and for coefficient_2
				"movss	%10, %%xmm3\n\t"
				"movss	%11, %%xmm1\n\t"
				"unpcklps	%%xmm1, %%xmm3\n\t"
				// Now %xmm3 = [ * | * | coefficient_2[3] |
				// coefficient_2[0] ]
				"unpcklps	%%xmm3, %%xmm4\n\t"
				// Now %xmm4 = [ coefficient_2[3] | coefficient_1[3] |
				// coefficient_2[0] | coefficient_1[0] ]
				// Same as above, with %xmm4 <-> %xmm5 and loading
				// elements 1, 4
				"movss	%12, %%xmm5\n\t"
				"movss	%13, %%xmm3\n\t"
				"unpcklps	%%xmm3, %%xmm5\n\t"
				"movss	%14, %%xmm3\n\t"
				"movss	%15, %%xmm1\n\t"
				"unpcklps	%%xmm1, %%xmm3\n\t"
				"unpcklps	%%xmm3, %%xmm5\n\t"
				// Zero %xmm3, %xmm1
				"xorps	%%xmm3, %%xmm3\n\t"
				"xorps	%%xmm1, %%xmm1"
				: "=m" (coefficient_1[0]), "=m" (coefficient_1[3]),
				  "=m" (coefficient_2[0]), "=m" (coefficient_2[3]),
				  "=m" (coefficient_1[1]), "=m" (coefficient_1[4]),
				  "=m" (coefficient_2[1]), "=m" (coefficient_2[4])
				: "m" (coefficient_1[0]), "m" (coefficient_1[3]),
				  "m" (coefficient_2[0]), "m" (coefficient_2[3]),
				  "m" (coefficient_1[1]), "m" (coefficient_1[4]),
				  "m" (coefficient_2[1]), "m" (coefficient_2[4])
				: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
				  "%xmm5");

			float *pinput = input;

			for(int i = n >> 1; i > 0; i--) {
				__asm__ __volatile__ (
					// Start with:
					// %xmm3 = [ * | * | y2[i - 1] | y1[i - 1] ]
					// %xmm1 = [ y2[i - 2] | y1[i - 2] | x[i - 1] |
					// x[i - 1] ]
					// Load x[i] into %xmm0
					"movss	%1, %%xmm0\n\t"
					// %xmm0 = [ * | * | * | x[i] ]
					"shufps	$0x40, %%xmm3, %%xmm0\n\t"
					// %xmm0 = [ y2[i - 1] | y1[i - 1] | x[i] | x[i] ]
					// %xmm1 = [ y2[i - 2] | y1[i - 2] | x[i - 1] |
					// x[i - 1] ]
					// Continue analogously to the coefficient[0] == 0
					// case.
					"movaps	%%xmm0, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"movaps	%%xmm1, %%xmm3\n\t"
					"mulps	%%xmm5, %%xmm3\n\t"
					"addps	%%xmm2, %%xmm3\n\t"
					"movhlps	%%xmm3, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm3"
					// %xmm3 = [ * | * | y2[i] | y1[i] ]
					: "=m" (*pinput)
					: "m" (*pinput)
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5");
				pinput++;
				__asm__ __volatile__ (
					// %xmm0 <-> %xmm1
					"movss	%1, %%xmm1\n\t"
					"shufps	$0x40, %%xmm3, %%xmm1\n\t"
					"movaps	%%xmm1, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"movaps	%%xmm0, %%xmm3\n\t"
					"mulps	%%xmm5, %%xmm3\n\t"
					"addps	%%xmm2, %%xmm3\n\t"
					"movhlps	%%xmm3, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm3"
					: "=m" (*pinput)
					: "m" (*pinput)
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5");
				pinput++;
			}
			pinput = input;

			float *poutput = output;

			for(int i = n >> 1; i > 0; i--) {
				__asm__ __volatile__ (
					// Start analogously as above.
					"movss	%2, %%xmm0\n\t"
					"shufps	$0x40, %%xmm3, %%xmm0\n\t"
					"movaps	%%xmm0, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"movaps	%%xmm1, %%xmm3\n\t"
					"mulps	%%xmm5, %%xmm3\n\t"
					"addps	%%xmm2, %%xmm3\n\t"
					"movhlps	%%xmm3, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm3\n\t"
					// %xmm3 = [ * | * | y2[i] | y1[i] ]
					"movaps	%%xmm3, %%xmm2\n\t"
#ifdef HAVE_SSE3
					"haddps	%%xmm2, %%xmm2\n\t"
#else // HAVE_SSE3
					"shufps	$0x55, %%xmm2, %%xmm2\n\t"
					"addss	%%xmm3, %%xmm2\n\t"
#endif // HAVE_SSE3
					"movss	%%xmm2, %1"
					: "=m" (*pinput), "=m" (*poutput)
					: "m" (*pinput), "m" (*poutput)
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5", "memory");
				pinput++;
				poutput++;
				__asm__ __volatile__ (
					// %xmm0 <-> %xmm1
					"movss	%2, %%xmm1\n\t"
					"shufps	$0x40, %%xmm3, %%xmm1\n\t"
					"movaps	%%xmm1, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"movaps	%%xmm0, %%xmm3\n\t"
					"mulps	%%xmm5, %%xmm3\n\t"
					"addps	%%xmm2, %%xmm3\n\t"
					"movhlps	%%xmm3, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm3\n\t"
					"movaps	%%xmm3, %%xmm2\n\t"
#ifdef HAVE_SSE3
					"haddps	%%xmm2, %%xmm2\n\t"
#else // HAVE_SSE3
					"shufps	$0x55, %%xmm2, %%xmm2\n\t"
					"addss	%%xmm3, %%xmm2\n\t"
#endif // HAVE_SSE3
					"movss	%%xmm2, %1"
					: "=m" (*pinput), "=m" (*poutput)
					: "m" (*pinput), "m" (*poutput)
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5");
				pinput++;
				poutput++;
			}
		}

#ifndef NVERIFY
		ID(float x_1_i_minus_2 = 0.0F);
		ID(float x_1_i_minus_1 = 0.0F);
		ID(float y_1_i_minus_2 = 0.0F);
		ID(float y_1_i_minus_1 = 0.0F);
		ID(float x_2_i_minus_2 = 0.0F);
		ID(float x_2_i_minus_1 = 0.0F);
		ID(float y_2_i_minus_2 = 0.0F);
		ID(float y_2_i_minus_1 = 0.0F);
		if(std::fpclassify(coefficient_1[0]) == FP_ZERO &&
		   std::fpclassify(coefficient_2[0]) == FP_ZERO) {
			for(int i = 0; i < n; i++) {
				ID(const float x_1_i = input[i]);
				ID(const float y_1_i =
				   coefficient_1[1] * x_1_i_minus_1 +
				   coefficient_1[2] * x_1_i_minus_2 -
				   coefficient_1[3] * y_1_i_minus_1 -
				   coefficient_1[4] * y_1_i_minus_2);
				IS(y_1_i_minus_2 = y_1_i_minus_1);
				IS(y_1_i_minus_1 = y_1_i);
				IS(x_1_i_minus_2 = x_1_i_minus_1);
				IS(x_1_i_minus_1 = x_1_i);
				ID(const float x_2_i = input[i]);
				ID(const float y_2_i =
				   coefficient_2[1] * x_2_i_minus_1 +
				   coefficient_2[2] * x_2_i_minus_2 -
				   coefficient_2[3] * y_2_i_minus_1 -
				   coefficient_2[4] * y_2_i_minus_2);
				IS(y_2_i_minus_2 = y_2_i_minus_1);
				IS(y_2_i_minus_1 = y_2_i);
				IS(x_2_i_minus_2 = x_2_i_minus_1);
				IS(x_2_i_minus_1 = x_2_i);
			}
			for(int i = 0; i < n; i++) {
				ID(const float x_1_i = input[i]);
				ID(const float y_1_i =
				   coefficient_1[1] * x_1_i_minus_1 +
				   coefficient_1[2] * x_1_i_minus_2 -
				   coefficient_1[3] * y_1_i_minus_1 -
				   coefficient_1[4] * y_1_i_minus_2);
				IS(y_1_i_minus_2 = y_1_i_minus_1);
				IS(y_1_i_minus_1 = y_1_i);
				IS(x_1_i_minus_2 = x_1_i_minus_1);
				IS(x_1_i_minus_1 = x_1_i);
				ID(const float x_2_i = input[i]);
				ID(const float y_2_i =
				   coefficient_2[1] * x_2_i_minus_1 +
				   coefficient_2[2] * x_2_i_minus_2 -
				   coefficient_2[3] * y_2_i_minus_1 -
				   coefficient_2[4] * y_2_i_minus_2);
				IS(y_2_i_minus_2 = y_2_i_minus_1);
				IS(y_2_i_minus_1 = y_2_i);
				IS(x_2_i_minus_2 = x_2_i_minus_1);
				IS(x_2_i_minus_1 = x_2_i);
				I(FEQ(output[i], y_1_i + y_2_i));
			}
		}
		else {
			I(std::fpclassify(coefficient_1[2]) == FP_ZERO);
			I(std::fpclassify(coefficient_2[2]) == FP_ZERO);
			for(int i = 0; i < n; i++) {
				ID(const float x_1_i = input[i]);
				ID(const float y_1_i =
				   coefficient_1[0] * x_1_i +
				   coefficient_1[1] * x_1_i_minus_1 -
				   coefficient_1[3] * y_1_i_minus_1 -
				   coefficient_1[4] * y_1_i_minus_2);
				IS(y_1_i_minus_2 = y_1_i_minus_1);
				IS(y_1_i_minus_1 = y_1_i);
				IS(x_1_i_minus_1 = x_1_i);
				ID(const float x_2_i = input[i]);
				ID(const float y_2_i =
				   coefficient_2[0] * x_2_i +
				   coefficient_2[1] * x_2_i_minus_1 -
				   coefficient_2[3] * y_2_i_minus_1 -
				   coefficient_2[4] * y_2_i_minus_2);
				IS(y_2_i_minus_2 = y_2_i_minus_1);
				IS(y_2_i_minus_1 = y_2_i);
				IS(x_2_i_minus_1 = x_2_i);
			}
			for(int i = 0; i < n; i++) {
				ID(const float x_1_i = input[i]);
				ID(const float y_1_i =
				   coefficient_1[0] * x_1_i +
				   coefficient_1[1] * x_1_i_minus_1 -
				   coefficient_1[3] * y_1_i_minus_1 -
				   coefficient_1[4] * y_1_i_minus_2);
				IS(y_1_i_minus_2 = y_1_i_minus_1);
				IS(y_1_i_minus_1 = y_1_i);
				IS(x_1_i_minus_1 = x_1_i);
				ID(const float x_2_i = input[i]);
				ID(const float y_2_i =
				   coefficient_2[0] * x_2_i +
				   coefficient_2[1] * x_2_i_minus_1 -
				   coefficient_2[3] * y_2_i_minus_1 -
				   coefficient_2[4] * y_2_i_minus_2);
				IS(y_2_i_minus_2 = y_2_i_minus_1);
				IS(y_2_i_minus_1 = y_2_i);
				IS(x_2_i_minus_1 = x_2_i);
				I(FEQ(output[i], y_1_i + y_2_i));
			}
		}
#endif // NVERIFY
	}

	void reconstruction_filtering_iir_t::
	accumulate_backward_aos_2(float output[], float input[],
							  const int n, float *coefficient_1,
							  float *coefficient_2)
	{
		I(A(int i = 0, i < n, i++, F(input[i])));
		I(A(int i = 0, i < 5, i++, F(coefficient_1[i])));
		I(A(int i = 0, i < 5, i++, F(coefficient_2[i])));
#ifndef NVERIFY
		ID(float original[n]);
		IS(memcpy(original, output, n * sizeof(float)));
#endif // NVERIFY

		const int n_1 = n - 1;

		if(std::fpclassify(coefficient_1[0]) == FP_ZERO &&
		   std::fpclassify(coefficient_2[0]) == FP_ZERO) {
			__asm__ __volatile__ (
				"movss	%8, %%xmm4\n\t"
				"movss	%9, %%xmm3\n\t"
				"unpcklps	%%xmm3, %%xmm4\n\t"
				"movss	%10, %%xmm3\n\t"
				"movss	%11, %%xmm1\n\t"
				"unpcklps	%%xmm1, %%xmm3\n\t"
				"unpcklps	%%xmm3, %%xmm4\n\t"
				"movss	%12, %%xmm5\n\t"
				"movss	%13, %%xmm3\n\t"
				"unpcklps	%%xmm3, %%xmm5\n\t"
				"movss	%14, %%xmm3\n\t"
				"movss	%15, %%xmm1\n\t"
				"unpcklps	%%xmm1, %%xmm3\n\t"
				"unpcklps	%%xmm3, %%xmm5\n\t"
				"xorps	%%xmm0, %%xmm0\n\t"
				"xorps	%%xmm1, %%xmm1"
				: "=m" (coefficient_1[2]), "=m" (coefficient_1[4]),
				  "=m" (coefficient_2[2]), "=m" (coefficient_2[4]),
				  "=m" (coefficient_1[1]), "=m" (coefficient_1[3]),
				  "=m" (coefficient_2[1]), "=m" (coefficient_2[3])
				: "m" (coefficient_1[2]), "m" (coefficient_1[4]),
				  "m" (coefficient_2[2]), "m" (coefficient_2[4]),
				  "m" (coefficient_1[1]), "m" (coefficient_1[3]),
				  "m" (coefficient_2[1]), "m" (coefficient_2[3])
				: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
				  "%xmm5");

			float *pinput = input + n_1;

			for(int i = n >> 1; i > 0; i--) {
				__asm__ __volatile__ (
					"movaps	%%xmm0, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"movaps	%%xmm1, %%xmm3\n\t"
					"mulps	%%xmm5, %%xmm3\n\t"
					"addps	%%xmm2, %%xmm3\n\t"
					"movhlps	%%xmm3, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm3\n\t"
					"movss	%1, %%xmm0\n\t"
					"shufps	$0x40, %%xmm3, %%xmm0"
					: "=m" (*pinput)
					: "m" (*pinput)
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5");
				pinput--;
				__asm__ __volatile__ (
					// %xmm0 <-> %xmm1
					"movaps	%%xmm1, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"movaps	%%xmm0, %%xmm3\n\t"
					"mulps	%%xmm5, %%xmm3\n\t"
					"addps	%%xmm2, %%xmm3\n\t"
					"movhlps	%%xmm3, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm3\n\t"
					"movss	%1, %%xmm1\n\t"
					"shufps	$0x40, %%xmm3, %%xmm1"
					: "=m" (*pinput)
					: "m" (*pinput)
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5");
				pinput--;
			}
			pinput = input + n_1;

			float *poutput = output + n_1;

			for(int i = n >> 1; i > 0; i--) {
				__asm__ __volatile__ (
					"movaps	%%xmm0, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"movaps	%%xmm1, %%xmm3\n\t"
					"mulps	%%xmm5, %%xmm3\n\t"
					"addps	%%xmm2, %%xmm3\n\t"
					"movhlps	%%xmm3, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm3\n\t"
					"movaps	%%xmm3, %%xmm2\n\t"
#ifdef HAVE_SSE3
					"haddps	%%xmm2, %%xmm2\n\t"
#else // HAVE_SSE3
					"shufps	$0x55, %%xmm2, %%xmm2\n\t"
					"addss	%%xmm3, %%xmm2\n\t"
#endif // HAVE_SSE3
					// Accumulate y[i] += %xmm2
					"addss	%2, %%xmm2\n\t"
					"movss	%%xmm2, %0\n\t"
					"movss	%3, %%xmm0\n\t"
					"shufps	$0x40, %%xmm3, %%xmm0"
					: "=m" (*poutput), "=m" (*pinput)
					: "m" (*poutput), "m" (*pinput)
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5");
				pinput--;
				poutput--;
				__asm__ __volatile__ (
					// %xmm0 <-> %xmm1
					"movaps	%%xmm1, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"movaps	%%xmm0, %%xmm3\n\t"
					"mulps	%%xmm5, %%xmm3\n\t"
					"addps	%%xmm2, %%xmm3\n\t"
					"movhlps	%%xmm3, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm3\n\t"
					"movaps	%%xmm3, %%xmm2\n\t"
#ifdef HAVE_SSE3
					"haddps	%%xmm2, %%xmm2\n\t"
#else // HAVE_SSE3
					"shufps	$0x55, %%xmm2, %%xmm2\n\t"
					"addss	%%xmm3, %%xmm2\n\t"
#endif // HAVE_SSE3
					"addss	%2, %%xmm2\n\t"
					"movss	%%xmm2, %0\n\t"
					"movss	%3, %%xmm1\n\t"
					"shufps	$0x40, %%xmm3, %%xmm1"
					: "=m" (*poutput), "=m" (*pinput)
					: "m" (*poutput), "m" (*pinput)
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5");
				pinput--;
				poutput--;
			}
		}
		else {
			// assuming coefficient_1[2] == 0 && coefficient_2[2] == 0
			__asm__ __volatile__ (
				"movss	%8, %%xmm4\n\t"
				"movss	%9, %%xmm3\n\t"
				"unpcklps	%%xmm3, %%xmm4\n\t"
				"movss	%10, %%xmm3\n\t"
				"movss	%11, %%xmm1\n\t"
				"unpcklps	%%xmm1, %%xmm3\n\t"
				"unpcklps	%%xmm3, %%xmm4\n\t"
				"movss	%12, %%xmm5\n\t"
				"movss	%13, %%xmm3\n\t"
				"unpcklps	%%xmm3, %%xmm5\n\t"
				"movss	%14, %%xmm3\n\t"
				"movss	%15, %%xmm1\n\t"
				"unpcklps	%%xmm1, %%xmm3\n\t"
				"unpcklps	%%xmm3, %%xmm5\n\t"
				"xorps	%%xmm3, %%xmm3\n\t"
				"xorps	%%xmm1, %%xmm1"
				: "=m" (coefficient_1[0]), "=m" (coefficient_1[3]),
				  "=m" (coefficient_2[0]), "=m" (coefficient_2[3]),
				  "=m" (coefficient_1[1]), "=m" (coefficient_1[4]),
				  "=m" (coefficient_2[1]), "=m" (coefficient_2[4])
				: "m" (coefficient_1[0]), "m" (coefficient_1[3]),
				  "m" (coefficient_2[0]), "m" (coefficient_2[3]),
				  "m" (coefficient_1[1]), "m" (coefficient_1[4]),
				  "m" (coefficient_2[1]), "m" (coefficient_2[4])
				: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
				  "%xmm5");

			float *pinput = input + n_1;

			for(int i = n >> 1; i > 0; i--) {
				__asm__ __volatile__ (
					"movss	%1, %%xmm0\n\t"
					"shufps	$0x40, %%xmm3, %%xmm0\n\t"
					"movaps	%%xmm0, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"movaps	%%xmm1, %%xmm3\n\t"
					"mulps	%%xmm5, %%xmm3\n\t"
					"addps	%%xmm2, %%xmm3\n\t"
					"movhlps	%%xmm3, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm3"
					: "=m" (*pinput)
					: "m" (*pinput)
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5");
				pinput--;
				__asm__ __volatile__ (
					// %xmm0 <-> %xmm1
					"movss	%1, %%xmm1\n\t"
					"shufps	$0x40, %%xmm3, %%xmm1\n\t"
					"movaps	%%xmm1, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"movaps	%%xmm0, %%xmm3\n\t"
					"mulps	%%xmm5, %%xmm3\n\t"
					"addps	%%xmm2, %%xmm3\n\t"
					"movhlps	%%xmm3, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm3"
					: "=m" (*pinput)
					: "m" (*pinput)
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5");
				pinput--;
			}
			pinput = input + n_1;

			float *poutput = output + n_1;

			for(int i = n >> 1; i > 0; i--) {
				__asm__ __volatile__ (
					"movss	%2, %%xmm0\n\t"
					"shufps	$0x40, %%xmm3, %%xmm0\n\t"
					"movaps	%%xmm0, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"movaps	%%xmm1, %%xmm3\n\t"
					"mulps	%%xmm5, %%xmm3\n\t"
					"addps	%%xmm2, %%xmm3\n\t"
					"movhlps	%%xmm3, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm3\n\t"
					"movaps	%%xmm3, %%xmm2\n\t"
#ifdef HAVE_SSE3
					"haddps	%%xmm2, %%xmm2\n\t"
#else // HAVE_SSE3
					"shufps	$0x55, %%xmm2, %%xmm2\n\t"
					"addss	%%xmm3, %%xmm2\n\t"
#endif // HAVE_SSE3
					"addss	%3, %%xmm2\n\t"
					"movss	%%xmm2, %1"
					: "=m" (*pinput), "=m" (*poutput)
					: "m" (*pinput), "m" (*poutput)
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5");
				pinput--;
				poutput--;
				__asm__ __volatile__ (
					// %xmm0 <-> %xmm1
					"movss	%2, %%xmm1\n\t"
					"shufps	$0x40, %%xmm3, %%xmm1\n\t"
					"movaps	%%xmm1, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"movaps	%%xmm0, %%xmm3\n\t"
					"mulps	%%xmm5, %%xmm3\n\t"
					"addps	%%xmm2, %%xmm3\n\t"
					"movhlps	%%xmm3, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm3\n\t"
					"movaps	%%xmm3, %%xmm2\n\t"
#ifdef HAVE_SSE3
					"haddps	%%xmm2, %%xmm2\n\t"
#else // HAVE_SSE3
					"shufps	$0x55, %%xmm2, %%xmm2\n\t"
					"addss	%%xmm3, %%xmm2\n\t"
#endif // HAVE_SSE3
					"addss	%3, %%xmm2\n\t"
					"movss	%%xmm2, %1"
					: "=m" (*pinput), "=m" (*poutput)
					: "m" (*pinput), "m" (*poutput)
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5");
				pinput--;
				poutput--;
			}
		}

#ifndef NVERIFY
		ID(float x_1_i_plus_2 = 0.0F);
		ID(float x_1_i_plus_1 = 0.0F);
		ID(float y_1_i_plus_2 = 0.0F);
		ID(float y_1_i_plus_1 = 0.0F);
		ID(float x_2_i_plus_2 = 0.0F);
		ID(float x_2_i_plus_1 = 0.0F);
		ID(float y_2_i_plus_2 = 0.0F);
		ID(float y_2_i_plus_1 = 0.0F);
		if(std::fpclassify(coefficient_1[0]) == FP_ZERO &&
		   std::fpclassify(coefficient_2[0]) == FP_ZERO) {
			for(int i = n - 1; i >= 0; i--) {
				ID(const float x_1_i = input[i]);
				ID(const float y_1_i =
				   coefficient_1[1] * x_1_i_plus_1 +
				   coefficient_1[2] * x_1_i_plus_2 -
				   coefficient_1[3] * y_1_i_plus_1 -
				   coefficient_1[4] * y_1_i_plus_2);
				IS(y_1_i_plus_2 = y_1_i_plus_1);
				IS(y_1_i_plus_1 = y_1_i);
				IS(x_1_i_plus_2 = x_1_i_plus_1);
				IS(x_1_i_plus_1 = x_1_i);
				ID(const float x_2_i = input[i]);
				ID(const float y_2_i =
				   coefficient_2[1] * x_2_i_plus_1 +
				   coefficient_2[2] * x_2_i_plus_2 -
				   coefficient_2[3] * y_2_i_plus_1 -
				   coefficient_2[4] * y_2_i_plus_2);
				IS(y_2_i_plus_2 = y_2_i_plus_1);
				IS(y_2_i_plus_1 = y_2_i);
				IS(x_2_i_plus_2 = x_2_i_plus_1);
				IS(x_2_i_plus_1 = x_2_i);
			}
			for(int i = n - 1; i >= 0; i--) {
				ID(const float x_1_i = input[i]);
				ID(const float y_1_i =
				   coefficient_1[1] * x_1_i_plus_1 +
				   coefficient_1[2] * x_1_i_plus_2 -
				   coefficient_1[3] * y_1_i_plus_1 -
				   coefficient_1[4] * y_1_i_plus_2);
				IS(y_1_i_plus_2 = y_1_i_plus_1);
				IS(y_1_i_plus_1 = y_1_i);
				IS(x_1_i_plus_2 = x_1_i_plus_1);
				IS(x_1_i_plus_1 = x_1_i);
				ID(const float x_2_i = input[i]);
				ID(const float y_2_i =
				   coefficient_2[1] * x_2_i_plus_1 +
				   coefficient_2[2] * x_2_i_plus_2 -
				   coefficient_2[3] * y_2_i_plus_1 -
				   coefficient_2[4] * y_2_i_plus_2);
				IS(y_2_i_plus_2 = y_2_i_plus_1);
				IS(y_2_i_plus_1 = y_2_i);
				IS(x_2_i_plus_2 = x_2_i_plus_1);
				IS(x_2_i_plus_1 = x_2_i);
				I(FEQ(output[i], original[i] + y_1_i + y_2_i));
			}
		}
		else {
			I(std::fpclassify(coefficient_1[2]) == FP_ZERO);
			I(std::fpclassify(coefficient_2[2]) == FP_ZERO);
			for(int i = n - 1; i >= 0; i--) {
				ID(const float x_1_i = input[i]);
				ID(const float y_1_i =
				   coefficient_1[0] * x_1_i +
				   coefficient_1[1] * x_1_i_plus_1 -
				   coefficient_1[3] * y_1_i_plus_1 -
				   coefficient_1[4] * y_1_i_plus_2);
				IS(y_1_i_plus_2 = y_1_i_plus_1);
				IS(y_1_i_plus_1 = y_1_i);
				IS(x_1_i_plus_1 = x_1_i);
				ID(const float x_2_i = input[i]);
				ID(const float y_2_i =
				   coefficient_2[0] * x_2_i +
				   coefficient_2[1] * x_2_i_plus_1 -
				   coefficient_2[3] * y_2_i_plus_1 -
				   coefficient_2[4] * y_2_i_plus_2);
				IS(y_2_i_plus_2 = y_2_i_plus_1);
				IS(y_2_i_plus_1 = y_2_i);
				IS(x_2_i_plus_1 = x_2_i);
			}
			for(int i = n - 1; i >= 0; i--) {
				ID(const float x_1_i = input[i]);
				ID(const float y_1_i =
				   coefficient_1[0] * x_1_i +
				   coefficient_1[1] * x_1_i_plus_1 -
				   coefficient_1[3] * y_1_i_plus_1 -
				   coefficient_1[4] * y_1_i_plus_2);
				IS(y_1_i_plus_2 = y_1_i_plus_1);
				IS(y_1_i_plus_1 = y_1_i);
				IS(x_1_i_plus_1 = x_1_i);
				ID(const float x_2_i = input[i]);
				ID(const float y_2_i =
				   coefficient_2[0] * x_2_i +
				   coefficient_2[1] * x_2_i_plus_1 -
				   coefficient_2[3] * y_2_i_plus_1 -
				   coefficient_2[4] * y_2_i_plus_2);
				IS(y_2_i_plus_2 = y_2_i_plus_1);
				IS(y_2_i_plus_1 = y_2_i);
				IS(x_2_i_plus_1 = x_2_i);
				I(FEQ(output[i], original[i] + y_1_i + y_2_i));
			}
		}
#endif // NVERIFY
	}

	void reconstruction_filtering_iir_t::
	apply_forward_soa_4(float output[], float input[], const int n,
						float *coefficient, const int stride)
	{
		I((stride & 0x3) == 0);

		// Causal and anticausal Gaussian IIR have 4 nonzero
		// coefficients. Implemented carefully, the ring buffer and
		// coefficients will tightly fit into 8 XMM registers
		// available both on AMD64/Intel 64 and IA-32.

		// Load the coefficients into %xmm4 .. %xmm7
		if(std::fpclassify(coefficient[0]) == FP_ZERO)
			__asm__ __volatile__ (
				"movss	%2, %%xmm4\n\t"
				"shufps	$0x0, %%xmm4, %%xmm4\n\t"
				"movss	%3, %%xmm5\n\t"
				"shufps	$0x0, %%xmm5, %%xmm5"
				: "=m" (coefficient[1]), "=m" (coefficient[2])
				: "m" (coefficient[1]), "m" (coefficient[2])
				: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
				  "%xmm5", "%xmm6", "%xmm7");
		else	// assuming coefficient[2] == 0
			__asm__ __volatile__ (
				"movss	%2, %%xmm4\n\t"
				"shufps	$0x0, %%xmm4, %%xmm4\n\t"
				"movss	%3, %%xmm5\n\t"
				"shufps	$0x0, %%xmm5, %%xmm5"
				: "=m" (coefficient[0]), "=m" (coefficient[1])
				: "m" (coefficient[0]), "m" (coefficient[1])
				: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
				  "%xmm5", "%xmm6", "%xmm7");
		__asm__ __volatile__ (
			"movss	%2, %%xmm6\n\t"
			"shufps	$0x0, %%xmm6, %%xmm6\n\t"
			"movss	%3, %%xmm7\n\t"
			"shufps	$0x0, %%xmm7, %%xmm7"
			: "=m" (coefficient[3]), "=m" (coefficient[4])
			: "m" (coefficient[3]), "m" (coefficient[4])
			: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4", "%xmm5",
			  "%xmm6", "%xmm7");

		int i_stride = 0;

		__asm__ __volatile__ (
			// Zero %xmm0 .. %xmm3
			"xorps	%%xmm0, %%xmm0\n\t"
			"xorps	%%xmm1, %%xmm1\n\t"
			"xorps	%%xmm2, %%xmm2\n\t"
			"xorps	%%xmm3, %%xmm3"
			: :
			: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4", "%xmm5",
			  "%xmm6", "%xmm7");

		const int stride_2 = stride << 1;

		if(std::fpclassify(coefficient[0]) == FP_ZERO) {
			for(int i = (n >> 1); i > 0; i--) {
				__asm__ __volatile__ (
					// The common index (j : j + vl - 1 : 1) is
					// suppressed below.
					// Start with:
					// %xmm0 = x[i - 2]
					// %xmm1 = x[i - 1]
					// %xmm2 = y[i - 2]
					// %xmm3 = y[i - 1]
					// %xmm4 = coefficient[1]
					// %xmm5 = coefficient[2]
					// %xmm6 = coefficient[3]
					// %xmm7 = coefficient[4]
					// Evaluate %xmm0 = coefficient[2] * x[i - 2] -
					// coefficient[4] * y[i - 2] to free one register.
					"mulps	%%xmm5, %%xmm0\n\t"
					"mulps	%%xmm7, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm0\n\t"
					// Note: %xmm2 is now free
					// Now: %xmm0 = coefficient[2] * x[i - 2] -
					// coefficient[4] * y[i - 2]
					// Evaluate %xmm2 = coefficient[1] * x[i - 1]
					"movaps	%%xmm1, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					// Accumulate in %xmm0 (we have to do the %xmm0
					// <-> %xmm2 register exchage here, since the
					// other remaining subtraction is not symmetric)
					"addps	%%xmm0, %%xmm2\n\t"
					// Note: %xmm0 is now free
					// Now: %xmm2 = coefficient[1] * x[i - 1] +
					// coefficient[2] * x[i - 2] - coefficient[4] *
					// y[i - 2]
					// Evaluate %xmm0 = coefficient[3] * y[i - 1]
					"movaps	%%xmm3, %%xmm0\n\t"
					"mulps	%%xmm6, %%xmm0\n\t"
					// Now: %xmm2 = y[i]
					"subps	%%xmm0, %%xmm2\n\t"
					// Load x[i] into %xmm0
					"movaps	%2, %%xmm0\n\t"
					// Store %xmm2 = y[i]
					"movaps	%%xmm2, %1"
					// Now:
					// %xmm0 = x[i - 1]
					// %xmm1 = x[i - 2]
					// %xmm2 = y[i - 1]
					// %xmm3 = y[i - 2]
					: "=m" (input[i_stride]),
					  "=m" (output[i_stride])
					: "m" (input[i_stride]),
					  "m" (output[i_stride])
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5", "%xmm6", "%xmm7");
				__asm__ __volatile__ (
					// The same as above, with %xmm0 <-> %xmm1, %xmm2
					// <-> %xmm3
					"mulps	%%xmm5, %%xmm1\n\t"
					"mulps	%%xmm7, %%xmm3\n\t"
					"subps	%%xmm3, %%xmm1\n\t"
					"movaps	%%xmm0, %%xmm3\n\t"
					"mulps	%%xmm4, %%xmm3\n\t"
					"addps	%%xmm1, %%xmm3\n\t"
					"movaps	%%xmm2, %%xmm1\n\t"
					"mulps	%%xmm6, %%xmm1\n\t"
					"subps	%%xmm1, %%xmm3\n\t"
					"movaps	%2, %%xmm1\n\t"
					"movaps	%%xmm3, %1"
					: "=m" (input[i_stride + stride]),
					  "=m" (output[i_stride + stride])
					: "m" (input[i_stride + stride]),
					  "m" (output[i_stride + stride])
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5", "%xmm6", "%xmm7");
				i_stride += stride_2;
			}
			if((n & 0x1) != 0)
				__asm__ __volatile__ (
					// The first part of the unrolled loop above,
					// minus one redundant data movement
					"mulps	%%xmm5, %%xmm0\n\t"
					"mulps	%%xmm7, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm0\n\t"
					"movaps	%%xmm1, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"addps	%%xmm0, %%xmm2\n\t"
					"movaps	%%xmm3, %%xmm0\n\t"
					"mulps	%%xmm6, %%xmm0\n\t"
					"subps	%%xmm0, %%xmm2\n\t"
					"movaps	%%xmm2, %0"
					: "=m" (output[i_stride])
					: "m" (output[i_stride])
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5", "%xmm6", "%xmm7");
		}
		else {	// assuming coefficient[2] == 0
			for(int i = (n >> 1); i > 0; i--) {
				__asm__ __volatile__ (
					// The common index (j : j + vl - 1 : 1) is
					// suppressed below.
					// Start with:
					// %xmm0 = x[i - 2]
					// %xmm1 = x[i - 1]
					// %xmm2 = y[i - 2]
					// %xmm3 = y[i - 1]
					// %xmm4 = coefficient[0]
					// %xmm5 = coefficient[1]
					// %xmm6 = coefficient[3]
					// %xmm7 = coefficient[4]
					// Note: %xmm0 is now free
					// Evaluate %xmm2 = coefficient[4] * y[i - 2]
					"mulps	%%xmm7, %%xmm2\n\t"
					// Evaluate %xmm0 = coefficient[3] * y[i - 1]
					"movaps	%%xmm3, %%xmm0\n\t"
					"mulps	%%xmm6, %%xmm0\n\t"
					// Evaluate %xmm0 = coefficient[3] * y[i - 1] +
					// coefficient[4] * y[i - 2]
					"addps	%%xmm2, %%xmm0\n\t"
					// Note: %xmm2 is now free
					// We need to free another register here to
					// prevent a duplicate memory movement. The
					// candidate is %xmm1, which is going to become
					// x[i - 2] and therefore ignored in the
					// subsequent step.
					// Evaluate %xmm1 = coefficient[1] * x[i - 1] -
					// coefficient[3] * y[i - 1] - coefficient[4] *
					// y[i - 2]
					"mulps	%%xmm5, %%xmm1\n\t"
					"subps	%%xmm0, %%xmm1\n\t"
					// Note: %xmm0, %xmm2 is now free
					// Load %xmm0 = x[i]
					"movaps	%2, %%xmm0\n\t"
					"movaps	%%xmm0, %%xmm2\n\t"
					// Evaluate %xmm2 = coefficient[0] * x[i]
					"mulps	%%xmm4, %%xmm2\n\t"
					// Evaluate %xmm0 = coefficient[0] * x[i] +
					// coefficient[1] * x[i - 1] - coefficient[3] *
					// y[i - 1] - coefficient[4] * y[i - 2]
					"addps	%%xmm1, %%xmm2\n\t"
					"movaps	%%xmm2, %1"
					// Now:
					// %xmm0 = x[i - 1]
					// %xmm1 = x[i - 2]
					// %xmm2 = y[i - 1]
					// %xmm3 = y[i - 2]
					: "=m" (input[i_stride]),
					  "=m" (output[i_stride])
					: "m" (input[i_stride]),
					  "m" (output[i_stride])
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5", "%xmm6", "%xmm7");
				__asm__ __volatile__ (
					// The same as above, with %xmm0 <-> %xmm1, %xmm2
					// <-> %xmm3
					"mulps	%%xmm7, %%xmm3\n\t"
					"movaps	%%xmm2, %%xmm1\n\t"
					"mulps	%%xmm6, %%xmm1\n\t"
					"addps	%%xmm3, %%xmm1\n\t"
					"mulps	%%xmm5, %%xmm0\n\t"
					"subps	%%xmm1, %%xmm0\n\t"
					"movaps	%2, %%xmm1\n\t"
					"movaps	%%xmm1, %%xmm3\n\t"
					"mulps	%%xmm4, %%xmm3\n\t"
					"addps	%%xmm0, %%xmm3\n\t"
					"movaps	%%xmm3, %1"
					: "=m" (input[i_stride + stride]),
					  "=m" (output[i_stride + stride])
					: "m" (input[i_stride + stride]),
					  "m" (output[i_stride + stride])
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5", "%xmm6", "%xmm7");
				i_stride += stride_2;
			}
			if((n & 0x1) != 0)
				__asm__ __volatile__ (
					// The first part of the unrolled loop above
					"mulps	%%xmm7, %%xmm2\n\t"
					"movaps	%%xmm3, %%xmm0\n\t"
					"mulps	%%xmm6, %%xmm0\n\t"
					"addps	%%xmm2, %%xmm0\n\t"
					"mulps	%%xmm5, %%xmm1\n\t"
					"subps	%%xmm0, %%xmm1\n\t"
					"movaps	%2, %%xmm0\n\t"
					"movaps	%%xmm0, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"addps	%%xmm1, %%xmm2\n\t"
					"movaps	%%xmm2, %1"
					: "=m" (input[i_stride]),
					  "=m" (output[i_stride])
					: "m" (input[i_stride]),
					  "m" (output[i_stride])
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5", "%xmm6", "%xmm7");
		}
#ifndef NVERIFY
		const int vector_length = 4;

		if(std::fpclassify(coefficient[0]) == FP_ZERO)
			for(int j = 0; j < vector_length; j++) {
				ID(float x_i_minus_2 = 0.0F);
				ID(float x_i_minus_1 = 0.0F);
				ID(float y_i_minus_2 = 0.0F);
				ID(float y_i_minus_1 = 0.0F);
				for(int i = 0; i < n; i++) {
					ID(const float x_i = input[i * stride + j]);
					ID(const float y_i =
					   coefficient[1] * x_i_minus_1 +
					   coefficient[2] * x_i_minus_2 -
					   coefficient[3] * y_i_minus_1 -
					   coefficient[4] * y_i_minus_2);
					IS(y_i_minus_2 = y_i_minus_1);
					IS(y_i_minus_1 = y_i);
					IS(x_i_minus_2 = x_i_minus_1);
					IS(x_i_minus_1 = x_i);
					I(FEQ(output[i * stride + j], y_i));
				}
			}
		else {
			I(std::fpclassify(coefficient[2]) == FP_ZERO);
			for(int j = 0; j < vector_length; j++) {
				ID(float x_i_minus_1 = 0.0F);
				ID(float y_i_minus_2 = 0.0F);
				ID(float y_i_minus_1 = 0.0F);
				for(int i = 0; i < n; i++) {
					ID(const float x_i = input[i * stride + j]);
					ID(const float y_i =
					   coefficient[0] * x_i +
					   coefficient[1] * x_i_minus_1 -
					   coefficient[3] * y_i_minus_1 -
					   coefficient[4] * y_i_minus_2);
					IS(y_i_minus_2 = y_i_minus_1);
					IS(y_i_minus_1 = y_i);
					IS(x_i_minus_1 = x_i);
					I(FEQ(output[i * stride + j], y_i));
				}
			}
		}
#endif // NVERIFY
	}

	void reconstruction_filtering_iir_t::
	accumulate_forward_soa_4(float output[], float input[],
							 const int n, float *coefficient,
							 const int stride)
	{
		I((stride & 0x3) == 0);

#ifndef NVERIFY
		const int vector_length = 4;

		ID(float original[n * vector_length]);
		for(int i = 0; i < n; i++)
			IS(memcpy(original + i * vector_length,
					  output + i * stride,
					  vector_length * sizeof(float)));
#endif // NVERIFY

		if(std::fpclassify(coefficient[0]) == FP_ZERO)
			__asm__ __volatile__ (
				"movss	%2, %%xmm4\n\t"
				"shufps	$0x0, %%xmm4, %%xmm4\n\t"
				"movss	%3, %%xmm5\n\t"
				"shufps	$0x0, %%xmm5, %%xmm5"
				: "=m" (coefficient[1]), "=m" (coefficient[2])
				: "m" (coefficient[1]), "m" (coefficient[2])
				: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
				  "%xmm5", "%xmm6", "%xmm7");
		else	// assuming coefficient[2] == 0
			__asm__ __volatile__ (
				"movss	%2, %%xmm4\n\t"
				"shufps	$0x0, %%xmm4, %%xmm4\n\t"
				"movss	%3, %%xmm5\n\t"
				"shufps	$0x0, %%xmm5, %%xmm5"
				: "=m" (coefficient[0]), "=m" (coefficient[1])
				: "m" (coefficient[0]), "m" (coefficient[1])
				: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
				  "%xmm5", "%xmm6", "%xmm7");
		__asm__ __volatile__ (
			"movss	%2, %%xmm6\n\t"
			"shufps	$0x0, %%xmm6, %%xmm6\n\t"
			"movss	%3, %%xmm7\n\t"
			"shufps	$0x0, %%xmm7, %%xmm7"
			: "=m" (coefficient[3]), "=m" (coefficient[4])
			: "m" (coefficient[3]), "m" (coefficient[4])
			: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4", "%xmm5",
			  "%xmm6", "%xmm7");

		int i_stride = 0;

		__asm__ __volatile__ (
			"xorps	%%xmm0, %%xmm0\n\t"
			"xorps	%%xmm1, %%xmm1\n\t"
			"xorps	%%xmm2, %%xmm2\n\t"
			"xorps	%%xmm3, %%xmm3"
			: :
			: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4", "%xmm5",
			  "%xmm6", "%xmm7");

		const int stride_2 = stride << 1;

		if(std::fpclassify(coefficient[0]) == FP_ZERO) {
			for(int i = (n >> 1); i > 0; i--) {
				__asm__ __volatile__ (
					// Analogous to apply_forward_soa_4()
					"mulps	%%xmm5, %%xmm0\n\t"
					"mulps	%%xmm7, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm0\n\t"
					"movaps	%%xmm1, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"addps	%%xmm0, %%xmm2\n\t"
					"movaps	%%xmm3, %%xmm0\n\t"
					"mulps	%%xmm6, %%xmm0\n\t"
					"subps	%%xmm0, %%xmm2\n\t"
					// Note: %xmm0 is now free
					// Accumulate y[i] += %xmm2
					"movaps	%%xmm2, %%xmm0\n\t"
					"addps	%2, %%xmm0\n\t"
					"movaps	%%xmm0, %0\n\t"
					"movaps	%3, %%xmm0"
					: "=m" (output[i_stride]),
					  "=m" (input[i_stride])
					: "m" (output[i_stride]),
					  "m" (input[i_stride])
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5", "%xmm6", "%xmm7");
				__asm__ __volatile__ (
					// The same as above, with %xmm0 <-> %xmm1, %xmm2
					// <-> %xmm3
					"mulps	%%xmm5, %%xmm1\n\t"
					"mulps	%%xmm7, %%xmm3\n\t"
					"subps	%%xmm3, %%xmm1\n\t"
					"movaps	%%xmm0, %%xmm3\n\t"
					"mulps	%%xmm4, %%xmm3\n\t"
					"addps	%%xmm1, %%xmm3\n\t"
					"movaps	%%xmm2, %%xmm1\n\t"
					"mulps	%%xmm6, %%xmm1\n\t"
					"subps	%%xmm1, %%xmm3\n\t"
					"movaps	%%xmm3, %%xmm1\n\t"
					"addps	%2, %%xmm1\n\t"
					"movaps	%%xmm1, %0\n\t"
					"movaps	%3, %%xmm1"
					: "=m" (output[i_stride + stride]),
					  "=m" (input[i_stride + stride])
					: "m" (output[i_stride + stride]),
					  "m" (input[i_stride + stride])
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5", "%xmm6", "%xmm7");
				i_stride += stride_2;
			}
			if((n & 0x1) != 0)
				__asm__ __volatile__ (
					"mulps	%%xmm5, %%xmm0\n\t"
					"mulps	%%xmm7, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm0\n\t"
					"movaps	%%xmm1, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"addps	%%xmm0, %%xmm2\n\t"
					"movaps	%%xmm3, %%xmm0\n\t"
					"mulps	%%xmm6, %%xmm0\n\t"
					"subps	%%xmm0, %%xmm2\n\t"
					"addps	%1, %%xmm2\n\t"
					"movaps	%%xmm2, %0"
					: "=m" (output[i_stride])
					: "m" (output[i_stride])
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5", "%xmm6", "%xmm7");
		}
		else {	// assuming coefficient[2] == 0
			for(int i = (n >> 1); i > 0; i--) {
				__asm__ __volatile__ (
					// Analogous to apply_forward_soa_4()
					"mulps	%%xmm7, %%xmm2\n\t"
					"movaps	%%xmm3, %%xmm0\n\t"
					"mulps	%%xmm6, %%xmm0\n\t"
					"addps	%%xmm2, %%xmm0\n\t"
					"mulps	%%xmm5, %%xmm1\n\t"
					"subps	%%xmm0, %%xmm1\n\t"
					"movaps	%2, %%xmm0\n\t"
					"movaps	%%xmm0, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"addps	%%xmm1, %%xmm2\n\t"
					// Note: %xmm1 is now free
					// Accmulate y[i] += %xmm2
					"movaps	%%xmm2, %%xmm1\n\t"
					"addps	%3, %%xmm1\n\t"
					"movaps	%%xmm1, %1"
					: "=m" (input[i_stride]),
					  "=m" (output[i_stride])
					: "m" (input[i_stride]),
					  "m" (output[i_stride])
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5", "%xmm6", "%xmm7");
				__asm__ __volatile__ (
					// The same as above, with %xmm0 <-> %xmm1, %xmm2
					// <-> %xmm3
					"mulps	%%xmm7, %%xmm3\n\t"
					"movaps	%%xmm2, %%xmm1\n\t"
					"mulps	%%xmm6, %%xmm1\n\t"
					"addps	%%xmm3, %%xmm1\n\t"
					"mulps	%%xmm5, %%xmm0\n\t"
					"subps	%%xmm1, %%xmm0\n\t"
					"movaps	%2, %%xmm1\n\t"
					"movaps	%%xmm1, %%xmm3\n\t"
					"mulps	%%xmm4, %%xmm3\n\t"
					"addps	%%xmm0, %%xmm3\n\t"
					"movaps	%%xmm3, %%xmm0\n\t"
					"addps	%3, %%xmm0\n\t"
					"movaps	%%xmm0, %1"
					: "=m" (input[i_stride + stride]),
					  "=m" (output[i_stride + stride])
					: "m" (input[i_stride + stride]),
					  "m" (output[i_stride + stride])
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5", "%xmm6", "%xmm7");
				i_stride += stride_2;
			}
			if((n & 0x1) != 0)
				__asm__ __volatile__ (
					"mulps	%%xmm7, %%xmm2\n\t"
					"movaps	%%xmm3, %%xmm0\n\t"
					"mulps	%%xmm6, %%xmm0\n\t"
					"addps	%%xmm2, %%xmm0\n\t"
					"mulps	%%xmm5, %%xmm1\n\t"
					"subps	%%xmm0, %%xmm1\n\t"
					"movaps	%2, %%xmm0\n\t"
					"movaps	%%xmm0, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"addps	%%xmm1, %%xmm2\n\t"
					"addps	%3, %%xmm2\n\t"
					"movaps	%%xmm2, %1"
					: "=m" (input[i_stride]),
					  "=m" (output[i_stride])
					: "m" (input[i_stride]),
					  "m" (output[i_stride])
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5", "%xmm6", "%xmm7");
		}
#ifndef NVERIFY
		if(std::fpclassify(coefficient[0]) == FP_ZERO)
			for(int j = 0; j < vector_length; j++) {
				ID(float x_i_minus_2 = 0.0F);
				ID(float x_i_minus_1 = 0.0F);
				ID(float y_i_minus_2 = 0.0F);
				ID(float y_i_minus_1 = 0.0F);
				for(int i = 0; i < n; i++) {
					ID(const float x_i = input[i * stride + j]);
					ID(const float y_i =
					   coefficient[1] * x_i_minus_1 +
					   coefficient[2] * x_i_minus_2 -
					   coefficient[3] * y_i_minus_1 -
					   coefficient[4] * y_i_minus_2);
					IS(y_i_minus_2 = y_i_minus_1);
					IS(y_i_minus_1 = y_i);
					IS(x_i_minus_2 = x_i_minus_1);
					IS(x_i_minus_1 = x_i);
					I(FEQ(output[i * stride + j],
						  original[i * vector_length + j] + y_i));
				}
			}
		else {
			I(std::fpclassify(coefficient[2]) == FP_ZERO);
			for(int j = 0; j < vector_length; j++) {
				ID(float x_i_minus_1 = 0.0F);
				ID(float y_i_minus_2 = 0.0F);
				ID(float y_i_minus_1 = 0.0F);
				for(int i = 0; i < n; i++) {
					ID(const float x_i = input[i * stride + j]);
					ID(const float y_i =
					   coefficient[0] * x_i +
					   coefficient[1] * x_i_minus_1 -
					   coefficient[3] * y_i_minus_1 -
					   coefficient[4] * y_i_minus_2);
					IS(y_i_minus_2 = y_i_minus_1);
					IS(y_i_minus_1 = y_i);
					IS(x_i_minus_1 = x_i);
					I(FEQ(output[i * stride + j],
						  original[i * vector_length + j] + y_i));
				}
			}
		}
#endif // NVERIFY
	}

	void reconstruction_filtering_iir_t::
	accumulate_backward_soa_4(float output[], float input[],
							  const int n, float *coefficient,
							  const int stride)
	{
		I((stride & 0x3) == 0);

#ifndef NVERIFY
		const int vector_length = 4;

		ID(float original[n * vector_length]);
		for(int i = 0; i < n; i++)
			for(int j = 0; j < vector_length; j++)
				IS(original[i * vector_length + j] = output[i * stride + j]);
#endif // NVERIFY

		if(std::fpclassify(coefficient[0]) == FP_ZERO)
			__asm__ __volatile__ (
				"movss	%2, %%xmm4\n\t"
				"shufps	$0x0, %%xmm4, %%xmm4\n\t"
				"movss	%3, %%xmm5\n\t"
				"shufps	$0x0, %%xmm5, %%xmm5"
				: "=m" (coefficient[1]), "=m" (coefficient[2])
				: "m" (coefficient[1]), "m" (coefficient[2])
				: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
				  "%xmm5", "%xmm6", "%xmm7");
		else	// assuming coefficient[2] == 0
			__asm__ __volatile__ (
				"movss	%2, %%xmm4\n\t"
				"shufps	$0x0, %%xmm4, %%xmm4\n\t"
				"movss	%3, %%xmm5\n\t"
				"shufps	$0x0, %%xmm5, %%xmm5"
				: "=m" (coefficient[0]), "=m" (coefficient[1])
				: "m" (coefficient[0]), "m" (coefficient[1])
				: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
				  "%xmm5", "%xmm6", "%xmm7");
		__asm__ __volatile__ (
			"movss	%2, %%xmm6\n\t"
			"shufps	$0x0, %%xmm6, %%xmm6\n\t"
			"movss	%3, %%xmm7\n\t"
			"shufps	$0x0, %%xmm7, %%xmm7"
			: "=m" (coefficient[3]), "=m" (coefficient[4])
			: "m" (coefficient[3]), "m" (coefficient[4])
			: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4", "%xmm5",
			  "%xmm6", "%xmm7");

		int i_stride = (n - 1) * stride;

		__asm__ __volatile__ (
			"xorps	%%xmm0, %%xmm0\n\t"
			"xorps	%%xmm1, %%xmm1\n\t"
			"xorps	%%xmm2, %%xmm2\n\t"
			"xorps	%%xmm3, %%xmm3"
			: :
			: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4", "%xmm5",
			  "%xmm6", "%xmm7");

		const int stride_2 = stride << 1;

		if(std::fpclassify(coefficient[0]) == FP_ZERO) {
			for(int i = (n >> 1); i > 0; i--) {
				__asm__ __volatile__ (
					// Analogous to apply_forward_soa_4()
					"mulps	%%xmm5, %%xmm0\n\t"
					"mulps	%%xmm7, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm0\n\t"
					"movaps	%%xmm1, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"addps	%%xmm0, %%xmm2\n\t"
					"movaps	%%xmm3, %%xmm0\n\t"
					"mulps	%%xmm6, %%xmm0\n\t"
					"subps	%%xmm0, %%xmm2\n\t"
					// Note: %xmm0 is now free
					// Accumulate y[i] += %xmm2
					"movaps	%%xmm2, %%xmm0\n\t"
					"addps	%2, %%xmm0\n\t"
					"movaps	%%xmm0, %0\n\t"
					"movaps	%3, %%xmm0"
					: "=m" (output[i_stride]),
					  "=m" (input[i_stride])
					: "m" (output[i_stride]),
					  "m" (input[i_stride])
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5", "%xmm6", "%xmm7");
				__asm__ __volatile__ (
					// The same as above, with %xmm0 <-> %xmm1, %xmm2
					// <-> %xmm3
					"mulps	%%xmm5, %%xmm1\n\t"
					"mulps	%%xmm7, %%xmm3\n\t"
					"subps	%%xmm3, %%xmm1\n\t"
					"movaps	%%xmm0, %%xmm3\n\t"
					"mulps	%%xmm4, %%xmm3\n\t"
					"addps	%%xmm1, %%xmm3\n\t"
					"movaps	%%xmm2, %%xmm1\n\t"
					"mulps	%%xmm6, %%xmm1\n\t"
					"subps	%%xmm1, %%xmm3\n\t"
					"movaps	%%xmm3, %%xmm1\n\t"
					"addps	%2, %%xmm1\n\t"
					"movaps	%%xmm1, %0\n\t"
					"movaps	%3, %%xmm1"
					: "=m" (output[i_stride - stride]),
					  "=m" (input[i_stride - stride])
					: "m" (output[i_stride - stride]),
					  "m" (input[i_stride - stride])
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5", "%xmm6", "%xmm7");
				i_stride -= stride_2;
			}
			if((n & 0x1) != 0)
				__asm__ __volatile__ (
					"mulps	%%xmm5, %%xmm0\n\t"
					"mulps	%%xmm7, %%xmm2\n\t"
					"subps	%%xmm2, %%xmm0\n\t"
					"movaps	%%xmm1, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"addps	%%xmm0, %%xmm2\n\t"
					"movaps	%%xmm3, %%xmm0\n\t"
					"mulps	%%xmm6, %%xmm0\n\t"
					"subps	%%xmm0, %%xmm2\n\t"
					"addps	%1, %%xmm2\n\t"
					"movaps	%%xmm2, %0"
					: "=m" (output[i_stride])
					: "m" (output[i_stride])
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5", "%xmm6", "%xmm7");
		}
		else {	// assuming coefficient[2] == 0
			for(int i = (n >> 1); i > 0; i--) {
				__asm__ __volatile__ (
					// Analogous to apply_forward_soa_4()
					"mulps	%%xmm7, %%xmm2\n\t"
					"movaps	%%xmm3, %%xmm0\n\t"
					"mulps	%%xmm6, %%xmm0\n\t"
					"addps	%%xmm2, %%xmm0\n\t"
					"mulps	%%xmm5, %%xmm1\n\t"
					"subps	%%xmm0, %%xmm1\n\t"
					"movaps	%2, %%xmm0\n\t"
					"movaps	%%xmm0, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"addps	%%xmm1, %%xmm2\n\t"
					// Note: %xmm1 is now free
					// Accmulate y[i] += %xmm2
					"movaps	%%xmm2, %%xmm1\n\t"
					"addps	%3, %%xmm1\n\t"
					"movaps	%%xmm1, %1"
					: "=m" (input[i_stride]),
					  "=m" (output[i_stride])
					: "m" (input[i_stride]),
					  "m" (output[i_stride])
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5", "%xmm6", "%xmm7");
				__asm__ __volatile__ (
					// The same as above, with %xmm0 <-> %xmm1, %xmm2
					// <-> %xmm3
					"mulps	%%xmm7, %%xmm3\n\t"
					"movaps	%%xmm2, %%xmm1\n\t"
					"mulps	%%xmm6, %%xmm1\n\t"
					"addps	%%xmm3, %%xmm1\n\t"
					"mulps	%%xmm5, %%xmm0\n\t"
					"subps	%%xmm1, %%xmm0\n\t"
					"movaps	%2, %%xmm1\n\t"
					"movaps	%%xmm1, %%xmm3\n\t"
					"mulps	%%xmm4, %%xmm3\n\t"
					"addps	%%xmm0, %%xmm3\n\t"
					"movaps	%%xmm3, %%xmm0\n\t"
					"addps	%3, %%xmm0\n\t"
					"movaps	%%xmm0, %1"
					: "=m" (input[i_stride - stride]),
					  "=m" (output[i_stride - stride])
					: "m" (input[i_stride - stride]),
					  "m" (output[i_stride - stride])
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5", "%xmm6", "%xmm7");
				i_stride -= stride_2;
			}
			if((n & 0x1) != 0)
				__asm__ __volatile__ (
					"mulps	%%xmm7, %%xmm2\n\t"
					"movaps	%%xmm3, %%xmm0\n\t"
					"mulps	%%xmm6, %%xmm0\n\t"
					"addps	%%xmm2, %%xmm0\n\t"
					"mulps	%%xmm5, %%xmm1\n\t"
					"subps	%%xmm0, %%xmm1\n\t"
					"movaps	%2, %%xmm0\n\t"
					"movaps	%%xmm0, %%xmm2\n\t"
					"mulps	%%xmm4, %%xmm2\n\t"
					"addps	%%xmm1, %%xmm2\n\t"
					"addps	%3, %%xmm2\n\t"
					"movaps	%%xmm2, %1"
					: "=m" (input[i_stride]),
					  "=m" (output[i_stride])
					: "m" (input[i_stride]),
					  "m" (output[i_stride])
					: "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4",
					  "%xmm5", "%xmm6", "%xmm7");
		}
#ifndef NVERIFY
		if(std::fpclassify(coefficient[0]) == FP_ZERO)
			for(int j = 0; j < vector_length; j++) {
				ID(float x_i_plus_2 = 0.0F);
				ID(float x_i_plus_1 = 0.0F);
				ID(float y_i_plus_2 = 0.0F);
				ID(float y_i_plus_1 = 0.0F);
				for(int i = n - 1; i >= 0; i--) {
					ID(const float x_i = input[i * stride + j]);
					ID(const float y_i =
					   coefficient[1] * x_i_plus_1 +
					   coefficient[2] * x_i_plus_2 -
					   coefficient[3] * y_i_plus_1 -
					   coefficient[4] * y_i_plus_2);
					IS(y_i_plus_2 = y_i_plus_1);
					IS(y_i_plus_1 = y_i);
					IS(x_i_plus_2 = x_i_plus_1);
					IS(x_i_plus_1 = x_i);
					I(FEQ(output[i * stride + j],
						  original[i * vector_length + j] + y_i));
				}
			}
		else {
			I(std::fpclassify(coefficient[2]) == FP_ZERO);
			for(int j = 0; j < vector_length; j++) {
				ID(float x_i_plus_1 = 0.0F);
				ID(float y_i_plus_2 = 0.0F);
				ID(float y_i_plus_1 = 0.0F);
				for(int i = n - 1; i >= 0; i--) {
					ID(const float x_i = input[i * stride + j]);
					ID(const float y_i =
					   coefficient[0] * x_i +
					   coefficient[1] * x_i_plus_1 -
					   coefficient[3] * y_i_plus_1 -
					   coefficient[4] * y_i_plus_2);
					IS(y_i_plus_2 = y_i_plus_1);
					IS(y_i_plus_1 = y_i);
					IS(x_i_plus_1 = x_i);
					I(FEQ(output[i * stride + j],
						  original[i * vector_length + j] + y_i));
				}
			}
		}
#endif // NVERIFY
	}

#endif // HAVE_SSE
}
