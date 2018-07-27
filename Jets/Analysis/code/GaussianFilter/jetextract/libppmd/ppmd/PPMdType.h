// This file is part of PPMd project
// Written and distributed to public domain by Dmitry Shkarin 1997,
// 1999-2001, 2006
// Contents: compilation parameters and miscelaneous definitions
// Comments: system & compiler dependent file

#ifndef _PPMDTYPE_H_
#define _PPMDTYPE_H_
#define NDEBUG
#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <cassert>
#include <sstream>

#ifdef __LP64__
#define _64_NORMAL
#else
#define _32_NORMAL
#endif

#ifndef NDEBUG
BOOL TestCompilation();	// for testing our data types
#endif // NDEBUG

#define _USE_PREFETCHING	// it gives 2-6% speed gain

#if (GCC_VERSION < 3030)
#define __thread
#endif

namespace ppmd
{
    const uint32_t PPMdSignature = 0x84acaf8f;
    enum {
		PROG_VAR = 'J',
		MAX_O = 16	// maximum allowed model order
    };

    template<class T>
    inline T CLAMP(const T &X, const T &LoX, const T &HiX)
    {
		return (X >= LoX) ? ((X <= HiX) ? (X) : (HiX)) : (LoX);
    }
    template<class T>
    inline void SWAP(T &t1, T &t2)
    {
		T tmp = t1;
		t1 = t2;
		t2 = tmp;
    }

#if 0
	// PPMd module works with file streams via ...GETC/...PUTC macros
	// only
	typedef FILE _PPMD_FILE;
#define _PPMD_E_GETC(fp)	getc(fp)
#define _PPMD_E_PUTC(c, fp)	putc((c), fp)
#define _PPMD_D_GETC(fp)	getc(fp)
#define _PPMD_D_PUTC(c, fp)	putc((c), fp)
#else
	typedef std::stringstream _PPMD_FILE;

	template<typename S_t>
	inline int _PPMD_E_GETC(S_t pps)
	{
		return pps->get();
	}

	template<typename C_t, typename S_t>
	inline void _PPMD_E_PUTC(C_t c, S_t pps)
	{
		pps->put(c);
	}

	template<typename S_t>
	inline int _PPMD_D_GETC(S_t pps)
	{
		return pps->get();
	}

	template<typename C_t, typename S_t>
	inline void _PPMD_D_PUTC(C_t c, S_t pps)
	{
		pps->put(c);
	}
#endif

}

#endif // _PPMDTYPE_H_
