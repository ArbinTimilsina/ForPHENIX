// -*- mode: c++; -*-

#ifndef XJETREC_SOLVE_H_
#define XJETREC_SOLVE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <vector>
#include <cmath>
#include <iostream>
#include <jetbase/num.h>

/////////////////////////////////////////////////////////////////////
//
// RHIC PHENIX JET RECONSTRUCTION LIBRARY
//
// Numerics: Non-Linear Equation Solver

namespace jet {

	/////////////////////////////////////////////////////////////////

	// General support functions

	/**
	 * Evaluate a univariate function with parameter
	 *
	 * @param[in] x				variable
	 * @param[in] function		function to evaluate
	 * @param[in] coefficient	function coefficients
	 * @return					function value at x
	 */
	template<typename F_t, typename X_t, typename Coeff_t>
	inline F_t
	solve_1d_impl_evaluate(
		const X_t variable,
		F_t (*function)(const X_t,
						const std::vector<Coeff_t> &),
		const std::vector<Coeff_t> &coefficient)
	{
		return function(variable, coefficient);
	}

	/////////////////////////////////////////////////////////////////

	// Solve a univariate transcendental equation, Alefeld-Potra-Shi
	// method
	//
	// References:
	//
	// G. E. Alefeld, F. A. Potra, and Y. X. Shi, ACM Trans. Math.
	// Software 21, 327 (1995).
	//
	// J. Stoer, Numerische Mathematik, Vol. 1 (Springer, Berlin,
	// 1999).

	/**
	 * Obtain an approximate root of f(x) by modified Aitken-Neville
	 * interpolation at a, b, d, and e.
	 *
	 * @param[in] a		node and lower bound of the interval
	 * @param[in] b		node and upper bound of the interval
	 * @param[in] d		node in [a, b]
	 * @param[in] e		node in [a, b]
	 * @param[in] fa	f(a)
	 * @param[in] fb	f(b)
	 * @param[in] fe	f(d)
	 * @param[in] fe	f(e)
	 * @return			approximate root of f(x)
	 * @see				Alefeld et al. (1995), p. 333, Stoer (1999),
	 *					section 2.1.2
	 */
	double solve_1d_aps_impl_ipzero(
		double a, double b, double d, double e,
		double fa, double fb, double fd, double fe);

	/**
	 * Approximate the zero in (a, b) of the quadratic polynomial
	 * interpolating f(x) at a, b, and d using k Newton steps
	 *
	 * @param[in] a		lower bound of the interval
	 * @param[in] b		upper bound of the interval
	 * @param[in] d		node in [a, b]
	 * @param[in] fa	f(a)
	 * @param[in] fb	f(b)
	 * @param[in] fd	f(d)
	 * @param[in] k		number of Newton steps to take
	 * @return			approximate zero in (a,b) of the quadratic
	 *					polynomial
	 * @see				Alefeld et al. (1995), p. 330
	 */
	double solve_1d_aps_impl_newton_quadratic(
		const double a, const double b, const double d,
		const double fa, const double fb, const double fd,
		const int k);

	/**
	 * Determines the termination criterion, with tolerance = 2 * (2 *
	 * epsilon * |b|).
	 *
	 * @param[in] b			
	 * @param[in] epsilon	relative machine precision
	 * @return				termination criterion
	 * @see					Alefeld et al. (1995), p. 340
	 */
	inline double solve_1d_aps_impl_tolerance(
		const double b, const double epsilon)
	{
		return 4 * fabs(b) * epsilon;
	}

	/**
	 * Given current enclosing interval [a, b] and a number c in (a,
	 * b), if f(c) = 0 then sets the output a = c. Otherwise
	 * determines the new enclosing interval: [a, b] = [a, c] or [a,
	 * b] = [c, b]. Also updates the termination criterion
	 * corresponding to the new enclosing interval.
	 *
	 * @param[in] function			function to evaluate
	 * @param[in] coefficient		function coefficients
	 * @param[in, out] a			lower bound of the interval
	 * @param[in, out] b			upper bound of the interval
	 * @param[out] c				node, used to determine the new
	 *								enclosing interval
	 * @param[in, out] fa			f(a)
	 * @param[in, out] fb			f(b)
	 * @param[in, out] tolerance	termination criterion
	 * @param[out] d				node, if the new enclosing
	 *								interval is [a, c] then d = b,
	 *								otherwise d = a
	 * @param[out] fd f(d)
	 * @param[in] epsilon			relative machine precision
	 * @see							Alefeld et al. (1995), p. 330, 341
	 */
	void
	solve_1d_aps_impl_bracket(
		double (*function)(const double,
						   const std::vector<double_complex_t> &),
		const std::vector<double_complex_t> &coefficient,
		double &a, double &b, double &c, double &fa, double &fb,
		double &tolerance, double &d, double &fd,
		const double epsilon, unsigned int &nevaluation);

	/**
	 * Calculate the product of differences (fi - fj), i != j, with f1
	 * = f(a), f2 = f(b), f3 = f(d), f4 = f(e)
	 *
	 * @param fa	f(a)
	 * @param fb	f(b)
	 * @param fd	f(d)
	 * @param fe	f(e)
	 * @return		product of differences (fi - fj)
	 */
	inline double solve_1d_aps_impl_product_f(
		const double fa, const double fb,
		const double fd, const double fe)
	{
		return ((fa - fb) * (fa - fd) * (fa - fe) *
				(fb - fd) * (fb - fe) * (fd - fe));
	}

	/**
	 * Find either an exact solution or an approximate solution of the
	 * equation f(x) = 0 in the interval [a, b], using the
	 * Alefeld-Potra-Shi method.
	 *
	 * At the begining of each iteration, the current enclosing
	 * interval is recorded as [a0, b0]. The first iteration is simply
	 * a secant step. Starting with the second iteration, three steps
	 * are taken in each iteration. First two steps are either
	 * quadratic interpolation or cubic inverse interpolation. The
	 * third step is a double-size secant step. If the diameter of the
	 * enclosing interval obtained after those three steps is larger
	 * than 0.5 * (b0 - a0), then an additional bisection step will be
	 * taken.
	 *
	 * @param[in] function		function to evaluate
	 * @param[in] coefficient	function coefficients
	 * @param[in, out] a		lower bound of the interval
	 * @param[in, out] b		upper bound of the interval
	 * @param[in] epsilon		relative machine precision
	 * @return					solution of the equation
	 * @see						Alefeld et al. (1995), p. 335
	 */
	double
	solve_1d_aps_impl_rroot(
		double (*function)(const double,
						   const std::vector<double_complex_t> &),
		const std::vector<double_complex_t> &coefficient,
		double &a, double &b, const double epsilon,
		unsigned int &niteration, unsigned int &nevaluation,
		const bool verbose = false);

	/**
	 * Find either an exact solution or an approximate solution of the
	 * equation f(x) = 0 in the interval [a, b], using the
	 * Alefeld-Potra-Shi method (easy-to-use driver function for
	 * solve_1d_aps_impl_rroot()).
	 *
	 * @param[in] function		function to evaluate
	 * @param[in] coefficient	function coefficients
	 * @param[in, out] a		lower bound of the interval
	 * @param[in, out] b		upper bound of the interval
	 * @param[in] epsilon		relative machine precision
	 * @return					solution of the equation
	 */
	double
	solve_1d_aps(
		double (*function)(const double,
						   const std::vector<double_complex_t> &),
		const std::vector<double_complex_t> &coefficient,
		const double a, const double b, const double epsilon = 0);

}

#endif // XJETREC_SOLVE_H_
