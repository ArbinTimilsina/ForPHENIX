// -*- mode: c++; -*-

#ifndef JETBASE_NUM_H_
#define JETBASE_NUM_H_

// LinkDef files generated by rootcint do not include config.h
#if defined(__CINT__) || defined(G__DICTIONARY)
#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H
#endif // defined(__CINT__) || defined(G__DICTIONARY)

#include <vector>
#include <ostream>

#include <cstdio>
#include <cfloat>
#include <climits>
#include <cmath>
#include <stdint.h>
#ifndef __INTEL_COMPILER
#include <complex>
#endif // __INTEL_COMPILER

// If GNU NANA is used, C99 complex number is handled explicitly in
// dbc.h, as the "I" macro conflicts with NANA
#include <jetbase/dbc.h>

/////////////////////////////////////////////////////////////////////

#ifdef __CINT__
#ifndef NAN
// As of ROOT 5.20/00, rootcint fails to pick up NAN from cmath.
// PHENIX collaborators: Leave this definition as it is, this is the
// correct (and glibc compatible) definition!
namespace {
	const unsigned int __cint_nan_byte = 0x7fc00000UL;
}
#define NAN	(*reinterpret_cast<const float *>(&__cint_nan_byte))
#endif // NAN
#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif // M_PI
#endif // __CINT__

namespace {
#ifdef __INTEL_COMPILER
	typedef float _Complex float_complex_t;
	typedef double _Complex double_complex_t;
#else // __INTEL_COMPILER
#define __SIMULATE_C99_COMPLEX
	typedef std::complex<float> float_complex_t;
	typedef std::complex<double> double_complex_t;
#ifdef _Complex_I
#undef _Complex_I
#endif // _Complex_I
	const double_complex_t _Complex_I(0, 1);
	inline double creal(const double_complex_t z)
	{
		return z.real();
	}
	inline float crealf(const float_complex_t z)
	{
		return z.real();
	}
	inline double cimag(const double_complex_t z)
	{
		return z.imag();
	}
	inline float cimagf(const float_complex_t z)
	{
		return z.imag();
	}
	inline double cabs(const double_complex_t z)
	{
		return abs(z);
	}
	inline float cabsf(const float_complex_t z)
	{
		return abs(z);
	}
	inline double carg(const double_complex_t z)
	{
		return arg(z);
	}
	inline float cargf(const float_complex_t z)
	{
		return arg(z);
	}
	inline double_complex_t cexp(const double_complex_t z)
	{
		return exp(z);
	}
	inline float_complex_t cexpf(const float_complex_t z)
	{
		return exp(z);
	}
	inline double_complex_t operator+(double x, double_complex_t z)
	{
		return double_complex_t(x) + z;
	}
	inline double_complex_t operator-(double x, double_complex_t z)
	{
		return double_complex_t(x) - z;
	}
	inline double_complex_t operator*(double x, double_complex_t z)
	{
		return double_complex_t(x) * z;
	}
	inline double_complex_t operator/(double x, double_complex_t z)
	{
		return double_complex_t(x) / z;
	}
#endif // __INTEL_COMPILER
}

namespace jet {

	/////////////////////////////////////////////////////////////////

	//

	/////////////////////////////////////////////////////////////////

	// Template floating point bounds

#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 869)
#endif // __INTEL_COMPILER
	class machine_limit_t {
	public:
		static inline float spacing(float x)
		{
			return FLT_EPSILON;
		}

		static inline double spacing(double x)
		{
			return DBL_EPSILON;
		}

		static inline float min_magnitude(float x)
		{
			return FLT_MIN;
		}

		static inline double min_magnitude(double x)
		{
			return DBL_MIN;
		}

		static inline float max_magnitude(float x)
		{
			return FLT_MAX;
		}

		static inline double max_magnitude(double x)
		{
			return DBL_MAX;
		}

#ifndef __CINT__
		static inline float spacing(float_complex_t z)
		{
			return FLT_EPSILON;
		}

		static inline double spacing(double_complex_t z)
		{
			return DBL_EPSILON;
		}

		static inline float min_magnitude(float_complex_t z)
		{
			return FLT_MIN;
		}

		static inline double min_magnitude(double_complex_t z)
		{
			return DBL_MIN;
		}

		static inline float max_magnitude(float_complex_t z)
		{
			return FLT_MAX;
		}

		static inline double max_magnitude(double_complex_t z)
		{
			return DBL_MAX;
		}
#endif // __CINT__

		template<typename x_t>
		static inline double spacing(std::vector<x_t> x)
		{
			return spacing(x[0]);
		}

		/**
		 * Range for cartesian component without overflow.
		 */
		static inline float max_cartesian(float x)
		{
			// |x_i| <= sqrt((1 - 2^(-24)) * 2^128 / 3) =
			// 1.06502323392266711e+19 gives (3 * x_i^2) within the
			// numerical limit.
			return 1.06502323392266711e+19F;
		}

		/**
		 * Highest pseudorapidity representable with IEEE 754 single
		 * precision, critical to ensure the vector operations being
		 * completely floating point exception free.
		 */
		static inline float max_pseudorapidity(float x)
		{
			// Pseudorapidities with |eta| > atanh((1 - 2^(-24))
			// 2^128) = 89.415986232628298363 will produce an
			// overflow. Because of the way we calculate with
			// pseudorapidities, the limit |eta| <= atanh(1 - 2^(-24))
			// = 8.66433974209815538 usually found in naive
			// implementations is not present.
			return 89.415986232628298363F;
		}
	};
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif // __INTEL_COMPILER

	/////////////////////////////////////////////////////////////////

	// IEEE 754-2008 operations

	template<typename X_t>
	inline X_t minnum(const X_t x, const X_t y)
	{
		if(std::isnan(x)) {
			return y;
		}
		else if(std::isnan(y)) {
			return x;
		}
		else {
			return std::min(x, y);
		}
	}

	template<typename X_t>
	inline X_t maxnum(const X_t x, const X_t y)
	{
		if(std::isnan(x)) {
			return y;
		}
		else if(std::isnan(y)) {
			return x;
		}
		else {
			return std::max(x, y);
		}
	}

	template<typename X_t>
	inline X_t minnan(const X_t x, const X_t y)
	{
		if(std::isnan(x)) {
			return x;
		}
		else if(std::isnan(y)) {
			return y;
		}
		else {
			return std::min(x, y);
		}
	}

	template<typename X_t>
	inline X_t maxnan(const X_t x, const X_t y)
	{
		if(std::isnan(x)) {
			return x;
		}
		else if(std::isnan(y)) {
			return y;
		}
		else {
			return std::max(x, y);
		}
	}

	/////////////////////////////////////////////////////////////////

	// Unit in the last place calculation

	/**
	 * Returns the unit in the last place of a IEEE 754 single
	 * precision value x, ulp(x), with the definition
	 * |machine_number(x) - x| <= ulp(x) / 2 after Harrison and Muller
	 * and the extension ulp(NaN) := NaN by Kahan.
	 *
	 * @param[in] x IEEE 754 single precision value
	 * @return ulp(x)
	 */
	inline double ulp(const float x)
	{
		double ulp;
		int exponent;

		switch(std::fpclassify(x)) {
		case FP_INFINITE:
			// fall through
		case FP_NAN:
			ulp = x;
			break;
		case FP_ZERO:
			ulp = ldexp(1.0, DBL_MIN_EXP - DBL_MANT_DIG + 1);
			break;
		default:	// FP_NORMAL, FP_SUBNORMAL
			frexp(x, &exponent);
			exponent = std::max(DBL_MIN_EXP, exponent) -
				DBL_MANT_DIG + 1;
			ulp = exponent > DBL_MAX_EXP ? INFINITY :
				ldexp(1.0, exponent);
		}

		return ulp;
	}

	/**
	 * Returns the unit in the last place of a IEEE 754 single
	 * precision value x, ulp(x), with the definition
	 * |machine_number(x) - x| <= ulp(x) / 2 after Harrison and Muller
	 * and the extension ulp(NaN) := NaN by Kahan.
	 *
	 * @param[in] x IEEE 754 single precision value
	 * @return ulp(x)
	 */
	inline double ulpf(const float x)
	{
		double ulp;
		int exponent;

		switch(std::fpclassify(x)) {
		case FP_INFINITE:
			// fall through
		case FP_NAN:
			ulp = x;
			break;
		case FP_ZERO:
			ulp = ldexp(1.0, FLT_MIN_EXP - FLT_MANT_DIG + 1);
			break;
		default:	// FP_NORMAL, FP_SUBNORMAL
			frexpf(x, &exponent);
			exponent = std::max(FLT_MIN_EXP, exponent) -
				FLT_MANT_DIG + 1;
			ulp = exponent > FLT_MAX_EXP ? INFINITY :
				ldexp(1.0, exponent);
		}

		return ulp;
	}

	/////////////////////////////////////////////////////////////////

	// Conversion

	// Fast IEEE 754 floating point to fixpoint conversions

	// FIXME: Does it really work on big-endian machines?

	/**
	 * Convert the integer quantization floor(2^n * (x + offset) *
	 * scale) of a IEEE 754 single precision floating point value x,
	 * (-offset) <= x < (-offset + 1 / scale) by bitmasking.
	 *
	 * @param[in] x x
	 * @param[in] offset offset, equals the lower bound
	 * @param[in] scale scale, equals the inverse of the range
	 * @param[in] n number of leading bits to return
	 * @return floor(2^n * (x + offset) * scale)
	 */
	inline unsigned int ftofixp(const float x, const float offset,
								const float scale, const int n)
	{
		const float x1 = 1.0F + (x + offset) * scale;

		if(x1 < 1 || x1 >= 2)
			fprintf(stderr, "%s:%d: x1 = %g\n", __FILE__, __LINE__,
					x1);

		// Precondition
		I(x1 >= 1 && x1 < 2);

		const unsigned int mantissa_bits =
			*((unsigned int *)(&x1)) & ((1U << 23) - 1U);

		I(n <= 23);

		return mantissa_bits >> (23 - n);
	}

	/**
	 * Convert the integer quantization floor(2^n * (x + offset) *
	 * scale) of a IEEE 754 double precision floating point value x,
	 * (-offset) <= x < (-offset + 1 / scale) by bitmasking.
	 *
	 * @param[in] x x
	 * @param[in] offset offset, equals the lower bound
	 * @param[in] scale scale, equals the inverse of the range
	 * @param[in] n number of leading bits to return
	 * @return floor(2^n * (x + offset) * scale)
	 */
	inline unsigned long long int ftofixp(const double x,
										  const double offset,
										  const double scale, 
										  const int n)
	{
		const double x1 = 1 + (x + offset) * scale;

		I(x1 >= 1 && x1 < 2);

		const unsigned long long int mantissa_bits =
			*((unsigned long long int *)(&x1)) &
			((1ULL << 52) - 1ULL);

		I(n <= 52);

		return mantissa_bits >> (52 - n);
	}

	inline int approx_floor(const float x)
	{
		const float xa = x + (float)(3 << 21);

		return (*((int *)&xa) - 0x4ac00000) >> 1;
	}

	// Convert a floating point value into a string in standard
	// Mathematica ($MachinePrecision) notation

	/**
	 * Return the Mathematica $MachinePrecision C string
	 * representation of an IEEE 754 double precision number.
	 *
	 * @param[in] value to convert
	 * @return C string Mathematica $MachinePrecision representation
	 */
	char *mathematica_form(const double x);

	/////////////////////////////////////////////////////////////////

	// Fast integer operations

	/**
	 * Returns the absolut value of n, using arithmetic bitshift.
	 *
	 * @param n signed input value
	 * @return the absolute value of n
	 */
	inline int32_t arithmetic_abs(const int32_t n)
	{
		// No precondition

		const int32_t sign = (n >> 31);
		const int32_t retval = (n ^ sign) - sign;

		// Postcondition. Note the peculiar behaviour: abs(-INT_MIN)
		// == INT_MIN, yet INT_MIN == -INT_MIN still holds.
		IG(retval >= 0, n != INT_MIN);
		IG(retval == n, n >= 0);
		IG(retval == -n, n < 0);

		return retval;
	}

#ifndef __CINT__
	/**
	 * Swap the 32-bit signed integers n0, n1, without external
	 * storage and using pair-wise XOR.
	 *
	 * @param[in,out] n0 n0
	 * @param[in,out] n1 n1
	 */
	inline void arithmetic_swap(int32_t &n0, int32_t &n1)
	{
		ID(int32_t n0_pre = n0);
		ID(int32_t n1_pre = n1);

		// No precondition

		n0 ^= n1;
		n1 ^= n0;
		n0 ^= n1;

		// Postcondition
		I(n0 == n1_pre);
		I(n1 == n0_pre);
	}
#endif // __CINT__

	/////////////////////////////////////////////////////////////////

	// Sorting algorithms

	/**
	 * Sort the 3 value pairs (x0, y0), (x1, y1), (x2, y2) in the
	 * ascending order x0 <= x1 <= x2. By default std::swap() is used
	 * to exchange the entries.
	 *
	 * @param[in,out] x0
	 * @param[in,out] y0
	 * @param[in,out] x1
	 * @param[in,out] y1
	 * @param[in,out] x2
	 * @param[in,out] y2
	 * @swapx[in] optional, function used to swap the x entries
	 * @swapy[in] optional, function used to swap the y entries
	 */
	template<typename X_t, typename Y_t>
	inline void sort_3_2(X_t &x0, Y_t &y0, X_t &x1, Y_t &y1,
						 X_t &x2, Y_t &y2,
						 void (*swapx)(X_t &, X_t &) = std::swap,
						 void (*swapy)(Y_t &, Y_t &) = std::swap)
	{
		// Efficient 3 element sorting algorithm with an average case
		// complexity of Omega(5/2 (comparisons) + 7/6 (exchanges))

		// Notation: L = low, M = mid, H = high

		// No precondition

		if(x0 <= x1) {
			// x0 <= x1
			if(x1 <= x2) {
				// x0 <= x1 <= x2
				// nothing
			}
			else {
				// x0 <= x1, x2 < x1
				if(x0 <= x2) {
					// x0 <= x2 < x1
					swapx(x1, x2);	// {L, H, M} -> {L, M, H}
					swapy(y1, y2);
				}
				else {
					// x2 < x0 <= x1
					swapx(x1, x2);	// {M, H, L} -> {H, M, L}
					swapy(y1, y2);
					swapx(x0, x1);	// {H, M, L} -> {L, M, H}
					swapy(y0, y1);
				}
			}
		}
		else {
			// x1 < x0
			if(x0 <= x2) {
				// x1 <= x0 <= x2
				swapx(x0, x1);		// {M, L, H} -> {L, M, H}
				swapy(y0, y1);
			}
			else {
				// x1 < x0, x2 < x0
				if(x1 <= x2) {
					// x1 <= x2 < x0
					swapx(x0, x1);	// {H, L, M} -> {L, H, M}
					swapy(y0, y1);
					swapx(x1, x2);	// {L, H, M} -> {L, M, H}
					swapy(y1, y2);
				}
				else {
					// x2 < x1 < x0
					swapx(x0, x2);	// {H, M, L} -> {L, M, H}
					swapy(y0, y2);
				}
			}
		}

		// Postcondition
		I(x0 <= x1 && x1 <= x2);
	}

	/////////////////////////////////////////////////////////////////

	// Rational numbers

#ifndef __CINT__
	class rational_t {
	private:
		int _numerator;
		int _denominator;
		inline rational_t instance(const int numerator,
								   const int denominator) const
		{
			rational_t ret;

			ret._numerator = numerator;
			ret._denominator = denominator;

			return ret;
		}
		/**
		 * Return the greates common divisor gcd(u, v)
		 *
		 * @param[in] u u
		 * @param[in] v v
		 * @return gcd(u, v)
		 */
		int gcd(int u, int v) const;
		/**
		 * Return the rational approximation (u / v) of x
		 *
		 * @param[out] u u
		 * @param[out] v v
		 * @param[in] x x
		 * @param[in] n maximum length of continued fraction
		 * expansion, equals 96 by default
		 */
		void rational_approx(int &u, int &v, const float x,
							 const int n = 96) const;
		/**
		 * Return the rational approximation (u / v) of x
		 *
		 * @param[out] u u
		 * @param[out] v v
		 * @param[in] x x
		 * @param[in] n maximum length of continued fraction
		 * expansion, equals 96 by default
		 */
		void rational_approx(int &u, int &v, const double x,
							 const int n = 96) const;
	public:
		rational_t(void)
		{
		}
		rational_t(const int numerator, const int denominator = 1);
		rational_t(const double x)
		{
			rational_approx(_numerator, _denominator, x);
		}
		inline int numerator(void) const
		{
			return _numerator;
		}
		inline int &numerator(void)
		{
			return _numerator;
		}
		inline int denominator(void) const
		{
			return _denominator;
		}
		inline int &denominator(void)
		{
			return _denominator;
		}
		inline int sign(void) const
		{
			I(_denominator >= 0);

			return _numerator > 0 ? 1 : _numerator < 0 ? -1 : 0;
		}
		inline double numerical(void) const
		{
			return (double)_numerator / (double)_denominator;
		}

		rational_t operator+(const rational_t v) const;
		rational_t operator-(const rational_t v) const;
		rational_t operator*(const rational_t v) const;
		rational_t operator/(const rational_t v) const;
		rational_t operator+=(const rational_t v);
		rational_t operator-=(const rational_t v);
		rational_t operator*=(const rational_t v);
		rational_t operator/=(const rational_t v);
		friend rational_t operator-(const rational_t u);

		bool operator<(const rational_t v) const;
		bool operator>(const rational_t v) const;
		bool operator<=(const rational_t v) const;
		bool operator>=(const rational_t v) const;
		bool operator==(const rational_t v) const;
		bool operator!=(const rational_t v) const;

		bool operator<(const int v) const;
		bool operator>(const int v) const;
		bool operator<=(const int v) const;
		bool operator>=(const int v) const;
		bool operator==(const int v) const;
		bool operator!=(const int v) const;
		friend bool operator<(const int u, const rational_t v);
		friend bool operator>(const int u, const rational_t v);
		friend bool operator<=(const int u, const rational_t v);
		friend bool operator>=(const int u, const rational_t v);
		friend bool operator==(const int u, const rational_t v);
		friend bool operator!=(const int u, const rational_t v);

		inline bool operator<(const double v) const
		{
			return _numerator < _denominator * v;
		}

		inline bool operator>(const double v) const
		{
			return _numerator > _denominator * v;
		}

		inline bool operator<=(const double v) const
		{
			return _numerator <= _denominator * v;
		}

		inline bool operator>=(const double v) const
		{
			return _numerator >= _denominator * v;
		}
		friend bool operator<(const double u, const rational_t v);
		friend bool operator>(const double u, const rational_t v);
		friend bool operator<=(const double u, const rational_t v);
		friend bool operator>=(const double u, const rational_t v);
	};

	inline rational_t operator-(const rational_t u)
	{
		return u.instance(-u._numerator, u._denominator);
	}

	inline bool operator<(const double u, const rational_t v)
	{
		return u * v._denominator < v._numerator;
	}

	inline bool operator>(const double u, const rational_t v)
	{
		return u * v._denominator > v._numerator;
	}

	inline bool operator<=(const double u, const rational_t v)
	{
		return u * v._denominator <= v._numerator;
	}

	inline bool operator>=(const double u, const rational_t v)
	{
		return u * v._denominator >= v._numerator;
	}

	std::ostream &operator<<(std::ostream &os, const rational_t &u);
#endif // __CINT__

	/////////////////////////////////////////////////////////////////

	int log2(const unsigned int n);

	void invariant_unsigned_division(
		const unsigned int divisor, int &algorithm,
		unsigned int &multiplier, int &shift);

}

#endif // JETBASE_NUM_H_