// -*- mode: c++; -*-

#ifndef JETBASE_DBC_H_
#define JETBASE_DBC_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

//////////////////////////////////////////////////////////////////////

// Jet reconstruction numerics: design by contract

// Include C99 complex number header and dischard the "I" macro
#ifdef __INTEL_COMPILER
#include <complex>
#if defined(HAVE_GMM) && !(defined(__CINT__) || defined(G__DICTIONARY))
#include <gmm/gmm.h>
#endif // defined(HAVE_GMM) && !(defined(__CINT__) || defined(G__D...
#ifndef _Complex
#define _Complex __complex__
#endif // _Complex
#if !defined(_COMPLEX_H) && !(defined(__CINT__) || defined(G__DICTIONARY))
#include "/usr/include/complex.h"
#undef I
#endif // !defined(_COMPLEX_H) && !(defined(__CINT__) || defined(G...
#endif // __INTEL_COMPILER

#include <cstdio>
#include <cmath>

#ifdef HAVE_ROOT
// CERN ROOT headers incompatible with GNU NANA
#include <TMath.h>
#include <TLorentzVector.h>
#endif // HAVE_ROOT

#if defined(HAVE_NANA) && !(defined(__CINT__) || defined(G__DICTIONARY))
#include <nana.h>
#else // defined(HAVE_NANA) && !(defined(__CINT__) || defined(G__DICTIONARY))
// Fake GNU Nana for rootcint and when running without GNU Nana
// installed
#define I(expr)
#define IG(expr, guard)
#define ID(def)
#define IS(statement)
#define WITHOUT_NANA
#endif // defined(HAVE_NANA) && !(defined(__CINT__) || defined(G__DICTIONARY))

// Use GNU C compatible void cast

#ifndef __CINT__
#ifndef __ASSERT_VOID_CAST
#if defined __cplusplus && __GNUC_PREREQ (2,95)
# define __ASSERT_VOID_CAST static_cast<void>
#else // defined __cplusplus && __GNUC_PREREQ (2,95)
# define __ASSERT_VOID_CAST (void)
#endif // defined __cplusplus && __GNUC_PREREQ (2,95)
#endif // __ASSERT_VOID_CAST
#endif // __CINT__

// Extra macros

#define F(expr) \
    (::__assert_isfinite(expr))

#define FEQ(x, y) \
    (::__assert_testfeq((x), (y), ::__assert_ftol,	\
     __FILE__, __LINE__))

#define FLE(x, y) \
    (::__assert_testfle((x), (y), ::__assert_ftol,	\
     __FILE__, __LINE__))

#define FGE(x, y) \
    (::__assert_testfgsim((x), (y), ::__assert_ftol,\
     __FILE__, __LINE__))

#define FRANGE(x, lo, hi) \
    ((x) < (lo) ? FEQ((x), (lo)) :				\
     (x) > (hi) ? FEQ((x), (hi)) : 1)

#define ALIGNED(x) \
    (::__aligned((x), __FILE__, __LINE__))

#define __ASSERT_FTOL	2.44140625e-4	// pow(0.5, (23 + 1) / 2)

#define SETFTOL(t) \
    (__ASSERT_VOID_CAST (::__assert_ftol = (double)(t)))

#define MULFTOL(t) \
    (__ASSERT_VOID_CAST (::__assert_ftol *= (double)(t)))

#define RESETFTOL \
    (SETFTOL(__ASSERT_FTOL))

/////////////////////////////////////////////////////////////////////

// Intel C++ Compiler compatibility

#ifdef __INTEL_COMPILER
#undef A
#define A(i,c,n,a) /* ForAll */ \
    ({ \
	bool _A_result = true; \
	i; \
	while((c) && _A_result) { \
	    if(!(a)) { \
		_A_result = false; \
	    } \
	    n; \
	} \
	_A_result; \
    })
#endif // __INTEL_COMPILER

/////////////////////////////////////////////////////////////////////

// Utility functions

// FIXME: Use template limits

namespace {

    // Constants and global variables
    double __assert_ftol = __ASSERT_FTOL;

    template<typename _T>
    inline int __assert_isfinite(_T __f)
    {
		return std::isfinite(__f);
    }

    // Floating point equal (similar)
    inline bool __assert_testfeq(const double x, const double y,
								 const double tol, const char *__file,
								 const int __line)
    {
		const double d = x - y;
		const bool ret = ((d >= -tol && d <= tol) ||
						  (d >= -fabs(y) * tol &&
						   d <=  fabs(y) * tol));

		if(!ret) {
			fprintf(stderr, "%s:%d: %lg != %lg (tolerance: %g)\n",
					__file, __line, x, y, tol);
		}

		return ret;
    }

    // Floating point less equal (less similar)
    inline bool __assert_testfle(const double x, const double y,
								 const double tol,
								 const char *__file, const int __line)
    {
		const double d = x - y;
		const bool ret = (d <= tol || d <= fabs(y) * tol);

		if(!ret) {
			fprintf(stderr, "%s:%d: %lg > %lg (tolerance: %g)\n",
					__file, __line, x, y, tol);
		}

		return ret;
    }

    template<typename X_t>
    inline bool __aligned(X_t *x,
						  const char *__file, const int __line)
    {
		I(reinterpret_cast<unsigned long>(x) != 0);

#ifdef HAVE_SSE
		const bool ret =
			(reinterpret_cast<unsigned long>(x) & 0xfUL) == 0;
#else // HAVE_SSE
		const bool ret = true;
#endif // HAVE_SSE

		if(!ret) {
			fprintf(stderr, "%s:%d: unaligned %lu\n",
					__file, __line,
					reinterpret_cast<unsigned long>(x));
		}

		return ret;
    }
}

#endif // JETBASE_DBC_H_
