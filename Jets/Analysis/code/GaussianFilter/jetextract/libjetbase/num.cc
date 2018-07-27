#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H

#include <cstdlib>
#include <jetbase/num.h>

namespace jet {

	/////////////////////////////////////////////////////////////////

	// Compatibility

	// Fortran BLAS style machine constants

	inline long int i1mach(const long int i)
	{
		switch(i) {
		case 1:		return 5;	// standard input
		case 2:		return 6;	// standard output
		case 3:		return 7;	// standard punch
		case 4:		return 0;	// standard error
		case 5:		return 32;	// bits per integer
		case 6:		return sizeof(int);
		case 7:		return 2;	// base for integers
		case 8:		return 31;	// digits of integer base
		case 9:		return LONG_MAX;
		case 10:	return FLT_RADIX;
		case 11:	return FLT_MANT_DIG;
		case 12:	return FLT_MIN_EXP;
		case 13:	return FLT_MAX_EXP;
		case 14:	return DBL_MANT_DIG;
		case 15:	return DBL_MIN_EXP;
		case 16:	return DBL_MAX_EXP;
		}
		fprintf(stderr, "invalid argument: i1mach(%ld)\n", i);
		exit(1);

		return 0;	// some compilers demand return values
	}

	inline float r1mach(const long int i)
	{
		switch(i){
		case 1:		return FLT_MIN;
		case 2:		return FLT_MAX;
		case 3:		return FLT_EPSILON / FLT_RADIX;
		case 4:		return FLT_EPSILON;
		case 5:		return (float)log10((double)FLT_RADIX);
		}
		fprintf(stderr, "invalid argument: r1mach(%ld)\n", i);
		exit(1);

		return 0;	// else complaint of missing return value
	}

	inline double d1mach(const long int i)
	{
		switch(i) {
		case 1:		return DBL_MIN;
		case 2:		return DBL_MAX;
		case 3:		return DBL_EPSILON / FLT_RADIX;
		case 4:		return DBL_EPSILON;
		case 5:		return log10((double)FLT_RADIX);
		}
		fprintf(stderr, "invalid argument: d1mach(%ld)\n", i);
		exit(1);

		return 0;	// some compilers demand return values
	}

	/////////////////////////////////////////////////////////////////

	char *mathematica_form(const double x)
	{
		const int buflen = 32;
		char buf[buflen];
		char *retval = new char[buflen];

		switch(std::fpclassify(x)) {
		case FP_NAN:
			strcpy(retval, "Indeterminate");
			return retval;
		case FP_INFINITE:
			strcpy(retval, x >= 0 ? "Infinity" : "-Infinity");
			return retval;
		case FP_ZERO:
			strcpy(retval, "0.`");
			return retval;
		default:	// FP_NORMAL or FP_SUBNORMAL
			int i, exp = 0;

			snprintf(buf, buflen, "%.16e\n", x);
			for(i = 0; buf[i]; i++)
				if(buf[i] == 'e') {
					buf[i] = '\0';
					sscanf(buf + i + 1, "%d", &exp);
				}

			for(i--; i >= 0 && buf[i] == '0'; i--)
				buf[i] = '\0';

			if(exp > -6 && exp < 6) {
				// Range in which Mathematica uses decimal notation
				snprintf(buf, buflen, "%%.%df", 16 - exp);
				snprintf(retval, buflen, buf, x);
				for(i = strlen(retval) - 2;
					i >= 0 && retval[i] == '0'; i--)
					retval[i] = '\0';
				strcat(retval, "`");
			}
			else
				// Range in which Mathematica uses scientific notation
				snprintf(retval, buflen - 1, "%s`*^%d", buf, exp);

			// Verification

			ID(char check_buf[32]);
			ID(double x_verify);

			IS(memcpy(check_buf, retval, 32));
			IS(for(i = 0; check_buf[i]; i++)
				   if(check_buf[i] == '`' &&
					  check_buf[i + 1] == '*' &&
					  check_buf[i + 2] == '^') {
					   check_buf[i] = 'e';
					   for(int j = i + 1; retval[j + 1] != '\0'; j++)
						   check_buf[j] = retval[j + 2];
				   });
			IS(sscanf(check_buf, "%lf`", &x_verify));
			IG(FEQ(x_verify, x), std::isfinite(x));

			return retval;
		}
	}

	int rational_t::gcd(int u, int v) const
	{
		int ret;

		ID(const int initial_u = u);
		ID(const int initial_v = v);

		// Trivial and pathological cases (abs is not well defined for
		// either u or v == INT_MIN)
		if(u == 0 || u == INT_MIN)
			ret = std::abs(v);
		else if(v == 0 || v == INT_MIN)
			ret = std::abs(u);
		else {
			u = std::abs(u);
			v = std::abs(v);

			I(u >= 0);
			I(v >= 0);

			// Exponent
			int k = 0;

			while((u & 1) == 0 /* u is even */ &&
				  (v & 1) == 0 /* v is even */) {
				u >>= 1;	// u /= 2
				v >>= 1;	// v /= 2
				k++;
			}
			do {
				if((u & 1) == 0)	// u is even
					u >>= 1;		// u /= 2
				else if((v & 1) == 0)	// v is even
					v >>= 1;		// v /= 2
				else if(u >= v)		// u is odd, v is odd, u >= v
					u = (u - v) >> 1;
				else			// u is odd, v is odd, u < v
					v = (v - u) >> 1;
			} while(u > 0);

			ret = v << k;	// v * 2^k
		}

		// Postcondition. Note gcd(u, v) is undefined if either u or v
		// == INT_MIN.
		IG(ret <= std::abs(initial_u), initial_u != 0 &&
		   initial_u != INT_MIN);
		IG(ret <= std::abs(initial_v), initial_v != 0 &&
		   initial_v != INT_MIN);
		IG(initial_u % ret == 0, ret != 0 && initial_u != INT_MIN);
		IG(initial_v % ret == 0, ret != 0 && initial_v != INT_MIN);
		IG(ret == std::abs(v), initial_u == 0);
		IG(ret == std::abs(u), initial_v == 0);

		return ret;
	}

	void rational_t::rational_approx(int &u, int &v, const double x,
									 const int n) const
	{
		std::vector<int> ret;
		// Continued fraction expansion: initial r and uncertainty
		double r = x;
		double dr = 0.5 * DBL_EPSILON * fabs(x);
		// Wallis' method: initialization
		int u2 = 0;
		int v2 = 1;
		int u1 = 1;
		int v1 = 0;
		// Approximation error
		double e1 = DBL_MAX;

		for(int i = 0; i < n; i++) {
			// Continued fraction expansion: obtain coefficient
			const double a = floor(r);
			// Wallis' method: iteration step
			const int u0 = (int)(a * u1 + u2);
			const int v0 = (int)(a * v1 + v2);
			// Calculate new approximation error
			const double e0 = fabs((double)u0 / (double)v0 - x);

			// Detect integer overflow by an increase of approximation
			// error
			if(e0 >= e1)
				break;
			// Continued fraction expansion: save coefficient
			u = u0;
			v = v0;
			// Approximation stops exactly
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 1572)
#endif // __INTEL_COMPILER
			if(a == r)
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif // __INTEL_COMPILER
				break;
			// Continued fraction expansion: iteration step
			r = 1 / (r - a);
			// Update uncertainty and stop if it grows above 1
			dr *= fabs(r * r);
			if(dr >= 1)
				break;
			// Wallis' method: shift. Note that gcd(u0, v0) <= 1
			// follows from the property of the continuant polynomials
			// (e.g. Knuth, 1998, p. 357, 4.5.3-(8))
			u2 = u1;
			v2 = v1;
			u1 = u0;
			v1 = v0;
			// Approximation error: shift
			e1 = e0;
		}
	}

	void rational_t::rational_approx(int &u, int &v, const float x,
									 const int n) const
	{
		std::vector<int> ret;
		// Continued fraction expansion: initial r and uncertainty
		double r = x;
		double dr = 0.5 * FLT_EPSILON * fabs(x);
		// Wallis' method: initialization
		int u2 = 0;
		int v2 = 1;
		int u1 = 1;
		int v1 = 0;
		// Approximation error
		double e1 = FLT_MAX;

		for(int i = 0; i < n; i++) {
			// Continued fraction expansion: obtain coefficient
			const double a = floor(r);
			// Wallis' method: iteration step
			const int u0 = (int)(a * u1 + u2);
			const int v0 = (int)(a * v1 + v2);
			// Calculate new approximation error
			const double e0 = fabs((double)u0 / (double)v0 - x);

			// Detect integer overflow by an increase of approximation
			// error
			if(e0 >= e1)
				break;
			// Continued fraction expansion: save coefficient
			u = u0;
			v = v0;
			// Approximation stops exactly
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 1572)
#endif // __INTEL_COMPILER
			if(a == r)
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif // __INTEL_COMPILER
				break;
			// Continued fraction expansion: iteration step
			r = 1 / (r - a);
			// Update uncertainty and stop if it grows above 1
			dr *= fabs(r * r);
			if(dr >= 1)
				break;
			// Wallis' method: shift. Note that gcd(u0, v0) <= 1
			// follows from the property of the continuant polynomials
			// (e.g. Knuth, 1998, p. 357, 4.5.3-(8))
			u2 = u1;
			v2 = v1;
			u1 = u0;
			v1 = v0;
			// Approximation error: shift
			e1 = e0;
		}
	}

	rational_t::rational_t(const int numerator,
						   const int denominator)
	{
		if(denominator >= 0) {
			_numerator = numerator;
			_denominator = denominator;
		}
		else {
			_numerator = -numerator;
			_denominator = -denominator;
		}

		I(_denominator >= 0);

		if(_denominator > 1) {
			const int d = gcd(_numerator, _denominator);

			_numerator /= d;
			_denominator /= d;
		}

		I(FEQ(numerical(), (double)numerator /
			  (double)denominator));
	}

	rational_t rational_t::operator+(const rational_t v) const
	{
		const int d1 = gcd(_denominator, v._denominator);

		rational_t ret;

		if(d1 <= 1)
			ret = instance(_numerator * v._denominator +
						   v._numerator * _denominator,
						   _denominator * v._denominator);
		else {
			const int t = _denominator * (v._numerator / d1) +
				v._denominator * (_numerator / d1);
			const int d2 = gcd(t, d1);

			ret = instance(t / d2, (_denominator / d1) *
						   (v._denominator / d2));
		}

		I(FEQ(ret.numerical(), numerical() + v.numerical()));

		return ret;
	}

	rational_t rational_t::operator-(const rational_t v) const
	{
		const int d1 = gcd(_denominator, v._denominator);

		rational_t ret;

		if(d1 <= 1)
			ret = instance(_numerator * v._denominator -
						   v._numerator * _denominator,
						   _denominator * v._denominator);
		else {
			const int t = _denominator * (v._numerator / d1) -
				v._denominator * (_numerator / d1);
			const int d2 = gcd(t, d1);

			ret = instance(t / d2, (_denominator / d1) *
						   (v._denominator / d2));
		}

		I(FEQ(ret.numerical(), numerical() - v.numerical()));

		return ret;
	}

	rational_t rational_t::operator*(const rational_t v) const
	{
		const int d1 = gcd(_numerator, v._denominator);
		const int d2 = gcd(_denominator, v._numerator);

		const rational_t ret =
			instance((_numerator / d1) * (v._numerator / d2),
					 (_denominator / d2) *
					 (v._denominator / d1));

		I(FEQ(ret.numerical(), numerical() * v.numerical()));

		return ret;
	}

	rational_t rational_t::operator/(const rational_t v) const
	{
		const int d1 = gcd(_numerator, v._numerator);
		const int d2 = gcd(_denominator, v._denominator);

		const rational_t ret =
			instance((_numerator / d1) * (v._denominator / d2) *
					 v.sign(),
					 std::abs((_denominator / d2) *
							  (v._numerator / d1)));

		I(FEQ(ret.numerical(), numerical() / v.numerical()));

		return ret;
	}

	rational_t rational_t::operator+=(const rational_t v)
	{
		const int d1 = gcd(_denominator, v._denominator);

		ID(const rational_t u = *this);

		rational_t ret;

		if(d1 <= 1) {
			_numerator = _numerator * v._denominator +
				v._numerator * _denominator;
			_denominator *= v._denominator;
		}
		else {
			const int t = _denominator * (v._numerator / d1) +
				v._denominator * (_numerator / d1);
			const int d2 = gcd(t, d1);

			_numerator = t / d2;
			_denominator = (_denominator / d1) *
				(v._denominator / d2);
		}

		I(FEQ(this->numerical(), u.numerical() + v.numerical()));

		return *this;
	}
	rational_t rational_t::operator-=(const rational_t v)
	{
		const int d1 = gcd(_denominator, v._denominator);

		if(d1 <= 1) {
			_numerator = _numerator * v._denominator -
				v._numerator * _denominator;
			_denominator *= v._denominator;
		}
		else {
			const int t = _denominator * (v._numerator / d1) -
				v._denominator * (_numerator / d1);
			const int d2 = gcd(t, d1);

			_numerator = t / d2;
			_denominator = (_denominator / d1) *
				(v._denominator / d2);
		}

		I(FEQ(this->numerical(), numerical() - v.numerical()));

		return *this;
	}

	rational_t rational_t::operator*=(const rational_t v)
	{
		const int d1 = gcd(_numerator, v._denominator);
		const int d2 = gcd(_denominator, v._numerator);

		_numerator = (_numerator / d1) * (v._numerator / d2);
		_denominator = (_denominator / d2) *
			(v._denominator / d1);

		I(FEQ(this->numerical(), numerical() * v.numerical()));

		return *this;
	}

	rational_t rational_t::operator/=(const rational_t v)
	{
		const int d1 = gcd(_numerator, v._numerator);
		const int d2 = gcd(_denominator, v._denominator);

		_numerator = (_numerator / d1) * (v._denominator / d2) *
			v.sign();
		_denominator = std::abs((_denominator / d2) *
								(v._numerator / d1));

		I(FEQ(this->numerical(), numerical() / v.numerical()));

		return *this;
	}

	std::ostream &operator<<(std::ostream &os, const rational_t &u)
	{
		if(u.numerator() == 0)
			return os << "0";
		else if(u.numerator() > 0)
			return os << u.numerator() << " / " << u.denominator();
		else
			return os << "-(" << -u.numerator() << " / "
					  << u.denominator() << ')';
	}

	bool rational_t::operator<(const rational_t v) const
	{
		const int d1 = gcd(_numerator, v._numerator);
		const int d2 = gcd(_denominator, v._denominator);

		return (_numerator / d1) * (v._denominator / d2) <
			(_denominator / d2) * (v._numerator / d1);
	}

	bool rational_t::operator>(const rational_t v) const
	{
		const int d1 = gcd(_numerator, v._numerator);
		const int d2 = gcd(_denominator, v._denominator);

		return (_numerator / d1) * (v._denominator / d2) >
			(_denominator / d2) * (v._numerator / d1);
	}

	bool rational_t::operator<=(const rational_t v) const
	{
		const int d1 = gcd(_numerator, v._numerator);
		const int d2 = gcd(_denominator, v._denominator);

		return (_numerator / d1) * (v._denominator / d2) <=
			(_denominator / d2) * (v._numerator / d1);
	}

	bool rational_t::operator>=(const rational_t v) const
	{
		const int d1 = gcd(_numerator, v._numerator);
		const int d2 = gcd(_denominator, v._denominator);

		return (_numerator / d1) * (v._denominator / d2) >=
			(_denominator / d2) * (v._numerator / d1);
	}

	bool rational_t::operator==(const rational_t v) const
	{
		const int d1 = gcd(_numerator, v._numerator);
		const int d2 = gcd(_denominator, v._denominator);

		return (_numerator / d1) * (v._denominator / d2) ==
			(_denominator / d2) * (v._numerator / d1);
	}

	bool rational_t::operator!=(const rational_t v) const
	{
		const int d1 = gcd(_numerator, v._numerator);
		const int d2 = gcd(_denominator, v._denominator);

		return (_numerator / d1) * (v._denominator / d2) !=
			(_denominator / d2) * (v._numerator / d1);
	}

	bool rational_t::operator<(const int v) const
	{
		const int d1 = gcd(_numerator, v);

		return _numerator / d1 < _denominator * (v / d1);
	}

	bool rational_t::operator>(const int v) const
	{
		const int d1 = gcd(_numerator, v);

		return _numerator / d1 > _denominator * (v / d1);
	}

	bool rational_t::operator<=(const int v) const
	{
		const int d1 = gcd(_numerator, v);

		return _numerator / d1 <= _denominator * (v / d1);
	}

	bool rational_t::operator>=(const int v) const
	{
		const int d1 = gcd(_numerator, v);

		return _numerator / d1 >= _denominator * (v / d1);
	}

	bool rational_t::operator==(const int v) const
	{
		const int d1 = gcd(_numerator, v);

		return _numerator / d1 == _denominator * (v / d1);
	}

	bool rational_t::operator!=(const int v) const
	{
		const int d1 = gcd(_numerator, v);

		return _numerator / d1 != _denominator * (v / d1);
	}

	bool operator<(const int u, const rational_t v)
	{
		const int d1 = v.gcd(u, v._numerator);

		return (u / d1) * v._denominator < v._numerator / d1;
	}

	bool operator>(const int u, const rational_t v)
	{
		const int d1 = v.gcd(u, v._numerator);

		return (u / d1) * v._denominator > v._numerator / d1;
	}

	bool operator<=(const int u, const rational_t v)
	{
		const int d1 = v.gcd(u, v._numerator);

		return (u / d1) * v._denominator <= v._numerator / d1;
	}

	bool operator>=(const int u, const rational_t v)
	{
		const int d1 = v.gcd(u, v._numerator);

		return (u / d1) * v._denominator >= v._numerator / d1;
	}

	bool operator==(const int u, const rational_t v)
	{
		const int d1 = v.gcd(u, v._numerator);

		return (u / d1) * v._denominator == v._numerator / d1;
	}

	bool operator!=(const int u, const rational_t v)
	{
		const int d1 = v.gcd(u, v._numerator);

		return (u / d1) * v._denominator != v._numerator / d1;
	}

	/////////////////////////////////////////////////////////////////

	int log2(const unsigned int n)
	{
		unsigned int s = n;
		int log2_n = 0;

		s >>= 1;
		while(s != 0) {
			s >>= 1;
			log2_n++;
		}

		return log2_n;
	}

	void invariant_unsigned_division(
		const unsigned int divisor, int &algorithm,
		unsigned int &multiplier, int &shift)
	{
		I(divisor < (1UL << 31));

		int exponent = 0;
		unsigned int odd = divisor;

		while((odd & 1UL) == 0) {
			odd >>= 1;
			exponent++;
		}

		if(odd == 1) {
			if(exponent == 0) {
				algorithm = 1;
				multiplier = 0xffffffffUL;
			}
			else {
				algorithm = 0;
				multiplier = 1UL << (32 - exponent);
			}
			shift = 0;
		}
		else {
			shift = log2(odd) + 1;

			unsigned long long k = (1ULL << (32 + shift)) /
				(0xffffffffULL - 0xffffffffULL %
				 (unsigned long long)odd);
			unsigned long long multiplier_low =
				(1ULL << (32 + shift)) / odd;
			unsigned long long multiplier_high =
				((1ULL << (32 + shift)) + k) / odd;

			while(((multiplier_low >> 1) < (multiplier_high >> 1)) &&
				  (shift > 0)) {
				multiplier_low >>= 1;
				multiplier_high >>= 1;
				shift--;
			}
			if((multiplier_high >> 32) == 0) {
				algorithm = 0;
				multiplier = (unsigned int)multiplier_high;
			}
			else {
				algorithm = 1;
				shift = log2(odd);
				multiplier_low =
					(1ULL << (32 + shift)) / (unsigned long long)odd;

				const unsigned int r =
					(unsigned int)((1ULL << (32 + shift)) %
								   (unsigned long long)odd);

				multiplier = r < ((odd >> 1) + 1) ?
					(unsigned int)multiplier_low :
					(unsigned int)multiplier_low + 1;
			}
		}
		while((multiplier & 1) == 0) {
			multiplier >>= 1;
			shift--;
		}
		shift += exponent;
	}
}
