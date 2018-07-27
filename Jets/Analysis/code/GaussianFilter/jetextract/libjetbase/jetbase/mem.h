// -*- mode: c++; -*-

#ifndef JETBASE_MEM_H_
#define JETBASE_MEM_H_

#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif // HAVE_CONFIG_H

/////////////////////////////////////////////////////////////////////

// Jet reconstruction numerics: memory allocation and copy

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600
#endif // _XOPEN_SOURCE
#include <cstdlib>
#include <cstring>

#include <jetbase/dbc.h>

namespace jet {

	/////////////////////////////////////////////////////////////////

	// Vector SIMD aligned memory allocation

	// Allocate memory for a general array, which will be 16 byte
	// aligned if SIMD is used (to speed up memcpy).
	inline void *malloc_aligned(const size_t n)
	{
		char *p = NULL;

		// Unfortunately, GCC 3.4 does not support the portable IA-32
		// _mm_malloc until very new versions. POSIX aligned malloc
		// used here.
		//
		// See also: http://gcc.gnu.org/bugzilla/show_bug.cgi?id=16570
#ifdef HAVE_SSE
#ifndef WITHOUT_NANA
		int ret = posix_memalign((void **)&p, 16, n);

		I(ret == 0);
#else // WITHOUT_NANA
		posix_memalign((void **)&p, 16, n);
#endif // WITHOUT_NANA

		I(ALIGNED(p));
#else // HAVE_SSE
		p = (char *)malloc(n);
#endif // HAVE_SSE

		return p;
	}

	// Allocate memory for an array of n Lorentz vectors. If SIMD
	// based vector operation is used, the returning pointer will
	// automatically satisfy the memory alignment required.
	inline float *malloc_lvec(const int n)
	{
		float *p = (float *)malloc_aligned(n * 4 * sizeof(float));

		I(p != NULL);

		return p;
	}

}

#ifdef HAVE_SSE

namespace jet {

	/////////////////////////////////////////////////////////////////

	// SSE vector SIMD based memory copy

	/**
	 * Copy __n bytes from memory area __src to __dest, using SSE
	 * instructions when appropriate
	 *
	 * @param __dest pointer to the destination memory area
	 * @param __src pointer to the source memory area
	 * @param __n number of bytes to copy
	 * @see memcpy(void *, const void *, size_t)
	 */
	void *memcpy_sse(void *__dest, const void *__src, size_t __n);

}

// Use SSE based memcpy by default hereafter
#define	memcpy(to, from, n)	jet::memcpy_sse((to), (from), (n))

#endif // HAVE_SSE

#endif // JETBASE_MEM_H_
