#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H

#include <cstdlib>
#include <cstring>
#include <stdint.h>
#include <jetbase/dbc.h>

namespace jet {

#ifdef HAVE_SSE
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 593)
#endif // __INTEL_COMPILER
    void *memcpy_sse(void *__dest, const void *__src, size_t __n)
    {
		// For efficiency, it is important that __src is 16 byte
		// aligned.

		// const unsigned long xmm_register_size = 16;

		// Precondition: valid pointers
		I(__src != NULL);
		I(__dest != NULL);

		// Precondition: no memory range overlap
		I(__src != __dest || __n == 0);
		IG((unsigned long)__src + __n <= (unsigned long)__dest,
		   (unsigned long)__src < (unsigned long)__dest);
		IG((unsigned long)__src >= (unsigned long)__dest + __n,
		   (unsigned long)__src > (unsigned long)__dest);

		ID(const size_t n_backup = __n);

		// C++ ideosyncracy
		uint8_t *src  = (uint8_t *)__src;
		uint8_t *dest = (uint8_t *)__dest;

		// FIXME: This cache prefetch is bad. For small __n, it
		// prefetches too much and unneccessarily through data out of
		// the cache. One should rather determine the amount to
		// prefetch based on __n.

		// Prefetch 320 bytes non-temoral data into the cache
		__asm__ __volatile__ (
			"prefetchnta	(%0)\n\t"
			"prefetchnta	64(%0)\n\t"
			"prefetchnta	128(%0)\n\t"
			"prefetchnta	192(%0)\n\t"
			"prefetchnta	256(%0)"
			:
			: "r" (src));

		if(__n >= 64 /* 4 * xmm_register_size */) {
			register unsigned long d = (unsigned long)dest &
				15 /* (xmm_register_size - 1) */;

			if(d != 0) {
				__n -= (d = 16 /* xmm_register_size */ - d);

				register unsigned long trash;

				// String based memcpy

				// Byte copy from [ESS] to [EDS], increment ESS, EDS,
				// decrement ECX until ECX == 0.
				__asm__ __volatile__ (
					"rep movsb"
					: "=&D" (dest), "=&S" (src), "=&c" (trash)
					: "0" (dest), "1" (src), "2" (d)
					: "memory");
			}
			if(((unsigned long)src & 15 /* (xmm_register_size - 1) */) != 0)
				// unaligned __src
				for(unsigned long i = __n >> 6 /* __n / 64 */; i > 0;
					i--) {
					// Prefetch additional 64 bytes into the cache.
					// Move 64 bytes unaligned from src to %xmm0 -
					// %xmm3. Store 64 bytes from %xmm0 - %xmm3 to
					// dest (unaligned) with non-temporal hint.
					__asm__ __volatile__ (
						"prefetchnta	320(%1)\n\t"
						"movups	(%1), %%xmm0\n\t"
						"movntps	%%xmm0, (%0)\n\t"
						"movups	16(%1), %%xmm1\n\t"
						"movntps	%%xmm1, 16(%0)\n\t"
						"movups	32(%1), %%xmm2\n\t"
						"movntps	%%xmm2, 32(%0)\n\t"
						"movups	48(%1), %%xmm3\n\t"
						"movntps	%%xmm3, 48(%0)"
						: "=&r" (dest)
						: "r" (src), "0" (dest)
						: "memory");
					src	 += 64 /* 4 * xmm_register_size */;
					dest += 64 /* 4 * xmm_register_size */;
				}
			else
				// aligned __src
				for(unsigned long i = __n >> 6 /* __n / 64 */; i > 0;
					i--) {
					// Prefetch additional 64 bytes into the cache.
					// Move 64 bytes aligned from src to %xmm0 -
					// %xmm3. Store 64 bytes from %xmm0 - %xmm3 to
					// dest (unaligned) with non-temporal hint.
					__asm__ __volatile__ (
						"prefetchnta	320(%1)\n\t"
						"movaps	(%1), %%xmm0\n\t"
						"movntps	%%xmm0, (%0)\n\t"
						"movaps	16(%1), %%xmm1\n\t"
						"movntps	%%xmm1, 16(%0)\n\t"
						"movaps	32(%1), %%xmm2\n\t"
						"movntps	%%xmm2, 32(%0)\n\t"
						"movaps	48(%1), %%xmm3\n\t"
						"movntps	%%xmm3, 48(%0)"
						: "=&r" (dest)
						: "r" (src), "0" (dest)
						: "memory");
					src	 += 64 /* 4 * xmm_register_size */;
					dest += 64 /* 4 * xmm_register_size */;
				}
			__n &= 63;	// __n = __n mod (4 * xmm_register_size)
		}
		if(__n != 0)
			// Use glibc memcpy() to copy the rest
			memcpy(dest, src, __n);

		// Postcondition: specification for the argument
		I(A(unsigned int i = 0, i < n_backup, i++,
			((char *)__dest)[i] == ((char *)__src)[i]));
		// Specification for return value is automatically fulfilled
		// by the return line. No postcondition checks necessary.

		return __dest;
    }
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif // __INTEL_COMPILER
#endif // HAVE_SSE

}
