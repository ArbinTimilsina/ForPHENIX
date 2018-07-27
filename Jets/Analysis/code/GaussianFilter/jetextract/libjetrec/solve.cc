#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H

#include <jetrec/solve.h>

namespace jet {

	double solve_1d_aps_impl_ipzero(
		double a, double b, double d, double e,
		double fa, double fb, double fd, double fe)
	{
#ifdef HAVE_SSE2
		// Vectorized version of solve_1d_aps_impl_ipzero() for SSE2
		double c;

		__asm__ __volatile__ (
			/////////////////////////////////////////////////////////
			// Part 1: Calculate:
			//
			// q21 = (b - d) * fb / (fd - fb);
			// q31 = (a - b) * fa / (fb - fa);
			// d21 = (b - d) * fd / (fd - fb);
			// d31 = (a - b) * fb / (fb - fa);
			//
			// Load %xmm0 = [ a | b ]
			"movhpd	%0, %%xmm0\n\t"
			"movlpd	%1, %%xmm0\n\t"
			// Load %xmm1 = [ b | d ]
			"movlpd	%2, %%xmm1\n\t"
			"shufpd	$0x0, %%xmm0, %%xmm1\n\t"
			// %xmm0 -= %xmm0
			"subpd	%%xmm1, %%xmm0\n\t"
			// %xmm0 = [ a - b | b - d ] now
			// Load %xmm3 = [ fa | fb ]
			"movsd	%4, %%xmm7\n\t"
			"unpcklpd	%%xmm7, %%xmm2\n\t"
			"movlpd	%5, %%xmm2\n\t"
			// Load %xmm2 = [ fb | fd ]
			"movlpd	%6, %%xmm1\n\t"
			"shufpd	$0x0, %%xmm2, %%xmm1\n\t"
			// Save [ fb | fd ] for the second interpolation step
			"movapd	%%xmm1, %%xmm3\n\t"
			// %xmm1 -= %xmm2
			"subpd	%%xmm2, %%xmm1\n\t"
			// %xmm1 = [ fb - fa | fd - fb ] now
			// %xmm2 = [ fa | fb ] is already in place and is being
			// reused.
			// %xmm2 *= %xmm0
			"mulpd	%%xmm0, %%xmm2\n\t"
			// %xmm2 = [ (a - b) * fa | (b - d) * fb ] now
			// %xmm2 /= %xmm1
			"divpd	%%xmm1, %%xmm2\n\t"
			// %xmm2 = [ (a - b) * fa / (fb - fa) | (b - d) * fb / (fd
			// - fb) ] now
			// Same procedure again for [ fb | fd ]
			"mulpd	%%xmm0, %%xmm3\n\t"
			"divpd	%%xmm1, %%xmm3\n\t"
			/////////////////////////////////////////////////////////
			// Part 2: Calculate:
			//
			// q11 = (d - e) * fd / (fe - fd);
			// q22 = (d21 - q11) * fb / (fe - fb);
			// q32 = (d31 - q21) * fa / (fd - fa);
			// d32 = (d31 - q21) * fd / (fd - fa);
			//
			// %xmm2 = [ q31 | q21 ], %xmm3 = [ d31 | d21 ] now
			"movsd	%2, %%xmm0\n\t"
			"subsd	%3, %%xmm0\n\t"
			"movsd	%7, %%xmm1\n\t"
			// Save %xmm5 = [ * | fe ]
			"movsd	%%xmm1, %%xmm5\n\t"
			"movsd	%6, %%xmm4\n\t"
			"subsd	%%xmm4, %%xmm1\n\t"
			// %xmm1 = [ * | fe - fd ] now
			"mulsd	%%xmm4, %%xmm0\n\t"
			"divsd	%%xmm1, %%xmm0\n\t"
			// %xmm0 = [ * | (d - e) * fd / (fe - fd) ] now
			// Now clogged: %xmm0 (q11), %xmm2 (q21, q31), %xmm3 (d21,
			// d31), %xmm4 (fd), %xmm5 (fe), %xmm7 (fa)
			"movsd	%5, %%xmm6\n\t"
			// %xmm6 = [ * | fb ] now
			"subsd	%%xmm6, %%xmm5\n\t"
			// %xmm5 = [ * | fe - fb ] now
			"movsd	%%xmm3, %%xmm1\n\t"
			// %xmm1 = [ * | d21 ] now
			"subsd	%%xmm0, %%xmm1\n\t"
			// %xmm1 = [ * | d21 - q11 ] now
			"mulsd	%%xmm6, %%xmm1\n\t"
			// %xmm1 = [ * | (d21 - q11) * fb ] now
			"divsd	%%xmm5, %%xmm1\n\t"
			// %xmm1 = [ * | (d21 - q11) * fb / (fe - fb) ] now
			// Now clogged: %xmm1 (q22), %xmm2 (q21, q31), %xmm3 (d21,
			// d31), %xmm5 (fe), %xmm4 (fd), %xmm7 (fa)
			"movapd	%%xmm3, %%xmm0\n\t"
			// %xmm1 = [ * | q22 ], %xmm2 = [ q31 | q21 ], %xmm3 = [
			// d31 | d21 ] now
			"unpckhpd	%%xmm3, %%xmm3\n\t"
			// %xmm3 = [ * | d31 ] now
			"subsd	%%xmm2, %%xmm3\n\t"
			// %xmm3 = [ * | d31 - q21 ] now
			"movsd	%%xmm4, %%xmm6\n\t"
			"subsd	%%xmm7, %%xmm6\n\t"
			// %xmm6 = [ * | fd - fa ] now
			"unpcklpd	%%xmm4, %%xmm4\n\t"
			"movsd	%%xmm7, %%xmm4\n\t"
			// %xmm4 = [ fd | fa ] now
			"unpcklpd	%%xmm3, %%xmm3\n\t"
			// %xmm3 = [ d31 - q21 | d31 - q21 ] now
			"mulpd	%%xmm4, %%xmm3\n\t"
			// %xmm3 = [ (d31 - q21) * fd | (d31 - q21) * fa ] now
			"unpcklpd	%%xmm6, %%xmm6\n\t"
			// %xmm6 = [ fd - fa | fd - fa ] now
			"divpd		%%xmm6, %%xmm3\n\t"
			// %xmm3 = [ (d31 - q21) * fd / (fd - fa) | (d31 - q21) *
			// fa / (fd - fa) ] now
			// Now clogged: %xmm1 (q22), %xmm2 (q31), %xmm3 (q32,
			// d32), either %xmm4 or %xmm7 (fa)
			/////////////////////////////////////////////////////////
			// Part 3: Calculate
			//
			// q33 = ((d32 - q22) * fa) / (fe - fa);
			"movsd	%7, %%xmm0\n\t"
			// %xmm0 = [ * | fe ] now
			"subsd	%%xmm7, %%xmm0\n\t"
			// %xmm0 = [ * | fe - fa ] now
			// %xmm3 = [ d32 | q32 ] now
			"movaps	%%xmm3, %%xmm4\n\t"
			"unpckhpd	%%xmm4, %%xmm4\n\t"
			// %xmm4 = [ * | d32 ] now
			"subsd	%%xmm1, %%xmm4\n\t"
			// %xmm4 = [ * | d32 - q22 ] now
			"mulsd	%%xmm4, %%xmm7\n\t"
			"divsd	%%xmm0, %%xmm7\n\t"
			// %xmm7 = [ * | q33 ] now
			/////////////////////////////////////////////////////////
			// Part 4: Calculate:
			//
			// c = q31 + q32 + q33 + a;
			"unpckhpd	%%xmm2, %%xmm2\n\t"
			// %xmm2 = [ * | q31 ] now
			// %xmm3 = [ * | q32 ] now
			"addsd	%%xmm3, %%xmm2\n\t"
			// %xmm1 = [ * | q31 + q32 ] now
			"addsd	%%xmm7, %%xmm2\n\t"
			// %xmm1 = [ * | q31 + q32 + q33 ] now
			"addsd	%0, %%xmm2\n\t"
			// %xmm1 = [ * | q31 + q32 + q33 + a ] now
			"movsd	%%xmm2, %8"
			: "=m" (a), "=m" (b), "=m" (d), "=m" (e),
			  "=m" (fa), "=m" (fb), "=m" (fd), "=m" (fe),
			  "=m" (c)
			: "m" (a), "m" (b), "m" (d), "m" (e),
			  "m" (fa), "m" (fb), "m" (fd), "m" (fe),
			  "m" (c)
			: "%xmm0", "%xmm1", "%xmm2", "%xmm3",
			  "%xmm4", "%xmm5", "%xmm6", "%xmm7");
#else // HAVE_SSE2
		double q11;
		double q21;
		double q31;
		double d21;
		double d31;

		q11 = (d - e) * fd / (fe - fd);
		q21 = (b - d) * fb / (fd - fb);
		q31 = (a - b) * fa / (fb - fa);
		d21 = (b - d) * fd / (fd - fb);
		d31 = (a - b) * fb / (fb - fa);

		double q22;
		double q32;
		double d32;

		q22 = (d21 - q11) * fb / (fe - fb);
		q32 = (d31 - q21) * fa / (fd - fa);
		d32 = (d31 - q21) * fd / (fd - fa);

		double q33;

		q33 = (d32 - q22) * fa / (fe - fa);

		// Calculate the output c

		const double c = q31 + q32 + q33 + a;
#endif // HAVE_SSE2

		return c;
	}

	double solve_1d_aps_impl_newton_quadratic(
		const double a, const double b, const double d,
		const double fa, const double fb, const double fd,
		const int k)
	{
		// Initialization: Find the coefficients of the quadratic
		// polynomial
		const double a0 = fa;
		const double a1 = (fb - fa) / (b - a);
		const double a2 = ((fd - fb) / (d - b) - a1) / (d - a);
		bool error = false;
		double c;

		do {
			// Safeguard to avoid overflow.
			if(std::fpclassify(a2) == FP_ZERO || error) {
				c = a - a0 / a1;
				return c;
			}
			// Determine the starting point of newton steps
			if((a2 > 0 && fa > 0) || (a2 < 0 && fa < 0))
				c = a;
			else
				c = b;
			// Start the safeguarded newton steps.
			for(int i = 0; i < k; i++) {
				if(!error) {
					const double pc =
						a0 + (a1 + a2 * (c - b)) * (c - a);
					const double pdc = a1 + a2 * ((2 * c) - (a + b));
					if(std::fpclassify(pdc) == FP_ZERO)
						error = true;
					else
						c -= pc / pdc;
				}
			}
		} while(error);

		return c;
	}

	void
	solve_1d_aps_impl_bracket(
		double (*function)(const double,
						   const std::vector<double_complex_t> &),
		const std::vector<double_complex_t> &coefficient,
		double &a, double &b, double &c, double &fa, double &fb,
		double &tolerance, double &d, double &fd,
		const double epsilon, unsigned int &nevaluation)
	{
		const double lambda = 0.7;

		// Adjust c if (b - a) is very small or if c is very close to
		// a or b
		tolerance *= lambda;
		if((b - a) <= 2 * tolerance)
			c = a + 0.5 * (b - a);
		else if(c <= a + tolerance)
			c = a + tolerance;
		else if(c >= b - tolerance)
			c = b - tolerance;

		// Call solve_1d_impl_evaluate() to obtain f(c)

		const double fc =
			solve_1d_impl_evaluate(c, function, coefficient);
		nevaluation++;

		// If f(c) = 0, then set a = c and return. This will terminate
		// the function in solve_1d_impl_rroot() and give the exact
		// solution of the equation f(x) = 0.
		if(std::fpclassify(fc) == FP_ZERO) {
			a = c;
			fa = 0;
			d = 0;
			fd = 0;
			return;
		}
		// If f(c) is not zero, then determine the new enclosing
		// interval
		if((fa > 0 && fc < 0) || (fa < 0 && fc > 0)) {
			d = b;
			fd = fb;
			b = c;
			fb = fc;
		}
		else {
			d = a;
			fd = fa;
			a = c;
			fa = fc;
		}
		// Update the termination criterion according to the new
		// enclosing interval
		tolerance =
			solve_1d_aps_impl_tolerance(fabs(fb) <= fabs(fa) ? b : a,
										epsilon);
	}

	double
	solve_1d_aps_impl_rroot(
		double (*function)(const double,
						   const std::vector<double_complex_t> &),
		const std::vector<double_complex_t> &coefficient,
		double &a, double &b, const double epsilon,
		unsigned int &niteration, unsigned int &nevaluation,
		const bool verbose)
	{
		const double mu = 0.5;

		// Initialization
		//
		// Set the number of iteration as 0. Call
		// solve_1d_impl_evaluate() to obtain the initial function
		// values f(a) and f(b). Set dumb values for the variable e
		// and fe.
		niteration = 0;

		double fa =
			solve_1d_impl_evaluate(a, function, coefficient);
		double fb =
			solve_1d_impl_evaluate(b, function, coefficient);
		double e = 1e+5;
		double fe = 1e+5;

		// 2 function evaluations so far
		nevaluation = 2;

		double c;
		double d;
		double fd;

		for(;;) {
			// Iteration starts. The enclosing interval before
			// executing the iteration is recorded as [a0, b0].
			const double a0 = a;
			const double b0 = b;

			// Updates the number of iteration.
			niteration++;

			// Calculates the termination criterion. Stops the
			// procedure if the criterion is satisfied.
			double tolerance = solve_1d_aps_impl_tolerance(
				fabs(fb) <= fabs(fa) ? b : a, epsilon);

			if(verbose) {
				std::cerr << __FILE__ << ':' << __LINE__
						  << ": information: " << "(a, b) = ("
						  << a << ", " << b << ')' << std::endl;
			}
			if((b - a) <= tolerance) {
				if(verbose) {
					std::cerr << __FILE__ << ':' << __LINE__
							  << ": information: " << "(a, b) = ("
							  << a << ", " << b << ')' << std::endl;
					std::cerr << __FILE__ << ':' << __LINE__
							  << ": information: " << "f(a) = " << fa
							  << std::endl;
					std::cerr << __FILE__ << ':' << __LINE__
							  << ": information: " << "b - a = "
							  << b - a << ", tolerance = "
							  << tolerance << std::endl;
				}
				// Terminates the procedure and return the root a.
				return a;
			}

			// For the first iteration, secant step is taken.
			if(niteration == 1) {
				c = a - (fa / (fb - fa)) * (b - a);

				// Call solve_1d_aps_impl_bracket() to get a shrinked
				// enclosing interval as well as to update the
				// termination criterion. stop the procedure if the
				// criterion is satisfied or the exact solution is
				// obtained.
				solve_1d_aps_impl_bracket(
					function, coefficient, a, b, c, fa, fb,
					tolerance, d, fd, epsilon, nevaluation);
				if(std::fpclassify(fa) == FP_ZERO ||
				   (b - a) <= tolerance) {
					if(verbose) {
						std::cerr << __FILE__ << ':' << __LINE__
								  << ": information: "
								  << "(a, b) = ("
								  << a << ", " << b << ')'
								  << std::endl;
						std::cerr << __FILE__ << ':' << __LINE__
								  << ": information: " << "f(a) = "
								  << fa << std::endl;
						std::cerr << __FILE__ << ':' << __LINE__
								  << ": information: " << "b - a = "
								  << b - a << ", tolerance = "
								  << tolerance << std::endl;
					}
					return a;
				}
				continue;
			}
			// Starting with the second iteration, in the first two
			// steps, either quadratic interpolation is used by
			// calling solve_1d_aps_impl_newton_quadratic() or the
			// cubic inverse interpolation is used by calling
			// solve_1d_aps_impl_ipzero(). In the following, if the
			// variable prof is not equal to 0, then the four function
			// value variables fa, fb, fd, and fe are distinct, and
			// hence solve_1d_aps_impl_ipzero() will be called.
			double prof =
				solve_1d_aps_impl_product_f(fa, fb, fd, fe);

			if(niteration == 2 || std::fpclassify(prof) == FP_ZERO)
				c = solve_1d_aps_impl_newton_quadratic(
					a, b, d, fa, fb, fd, 2);
			else {
				c = solve_1d_aps_impl_ipzero(
					a, b, d, e, fa, fb, fd, fe);
				if((c - a) * (c - b) >= 0)
					c = solve_1d_aps_impl_newton_quadratic(
						a, b, d, fa, fb, fd, 2);
			}

			e = d;
			fe = fd;

			// Call solve_1d_aps_impl_bracket() to get a shrinked
			// enclosing interval as well as to update the termination
			// criterion. Stop the procedure if the criterion is
			// satisfied or the exact solution is obtained.
			solve_1d_aps_impl_bracket(
				function, coefficient, a, b, c, fa, fb,
				tolerance, d, fd, epsilon, nevaluation);
			if(std::fpclassify(fa) == FP_ZERO ||
			   (b - a) <= tolerance) {
				if(verbose) {
					std::cerr << __FILE__ << ':' << __LINE__
							  << ": information: " << "(a, b) = ("
							  << a << ", " << b << ')' << std::endl;
					std::cerr << __FILE__ << ':' << __LINE__
							  << ": information: " << "f(a) = " << fa
							  << std::endl;
					std::cerr << __FILE__ << ':' << __LINE__
							  << ": information: " << "b - a = "
							  << b - a << ", tolerance = "
							  << tolerance << std::endl;
				}
				return a;
			}

			prof = solve_1d_aps_impl_product_f(fa, fb, fd, fe);
			if(std::fpclassify(prof) == FP_ZERO) {
				c = solve_1d_aps_impl_newton_quadratic(
					a, b, d, fa, fb, fd, 3);
			}
			else {
				c = solve_1d_aps_impl_ipzero(a, b, d, e,
											 fa, fb, fd, fe);
				if((c - a) * (c - b) >= 0)
					c = solve_1d_aps_impl_newton_quadratic(
						a, b, d, fa, fb, fd, 3);
			}
			// Call solve_1d_aps_impl_bracket() to get a shrinked
			// enclosing interval as well as to update the termination
			// criterion. Stop the procedure if the criterion is
			// satisfied or the exact solution is obtained.
			solve_1d_aps_impl_bracket(function, coefficient,
									  a, b, c, fa, fb,
									  tolerance, d, fd, epsilon,
									  nevaluation);
			if(std::fpclassify(fa) == FP_ZERO ||
			   (b - a) <= tolerance) {
				if(verbose) {
					std::cerr << __FILE__ << ':' << __LINE__
							  << ": information: " << "(a, b) = ("
							  << a << ", " << b << ')' << std::endl;
					std::cerr << __FILE__ << ':' << __LINE__
							  << ": information: " << "f(a) = " << fa
							  << std::endl;
					std::cerr << __FILE__ << ':' << __LINE__
							  << ": information: " << "b - a = "
							  << b - a << ", tolerance = "
							  << tolerance << std::endl;
				}
				return a;
			}

			e = d;
			fe = fd;

			double u;
			double fu;

			// Takes the double-size secant step.
			if(fabs(fa) < fabs(fb)) {
				u = a;
				fu = fa;
			}
			else {
				u = b;
				fu = fb;
			}
			c = u - 2 * (fu / (fb - fa)) * (b - a);
			if(fabs(c - u) > 0.5 * (b - a)) {
				c = a + 0.5 * (b - a);
			}
			// Call solve_1d_aps_impl_bracket() to get a shrinked
			// enclosing interval as well as to update the termination
			// criterion. Stop the procedure if the criterion is
			// satisfied or the exact solution is obtained.
			solve_1d_aps_impl_bracket(
				function, coefficient, a, b, c, fa, fb,
				tolerance, d, fd, epsilon, nevaluation);
			if(std::fpclassify(fa) == FP_ZERO ||
			   (b - a) <= tolerance) {
				if(verbose) {
					std::cerr << __FILE__ << ':' << __LINE__
							  << ": information: " << "(a, b) = ("
							  << a << ", " << b << ')' << std::endl;
					std::cerr << __FILE__ << ':' << __LINE__
							  << ": information: " << "f(a) = " << fa
							  << std::endl;
					std::cerr << __FILE__ << ':' << __LINE__
							  << ": information: " << "b - a = "
							  << b - a << ", tolerance = "
							  << tolerance << std::endl;
				}
				return a;
			}

			// Determines whether an additional bisection step is
			// needed. And takes it if necessary.
			if(b - a >= mu * (b0 - a0)) {
				e = d;
				fe = fd;

				// Call solve_1d_aps_impl_bracket() to get a shrinked
				// enclosing interval as well as to update the
				// termination criterion. Stop the procedure if the
				// criterion is satisfied or the exact solution is
				// obtained.
				double c_additional = a + 0.5 * (b - a);

				solve_1d_aps_impl_bracket(
					function, coefficient, a, b, c_additional, fa, fb,
					tolerance, d, fd, epsilon, nevaluation);
				if(std::fpclassify(fa) == FP_ZERO ||
				   (b - a) <= tolerance) {
					if(verbose) {
						std::cerr << __FILE__ << ':' << __LINE__
								  << ": information: "
								  << "(a, b) = ("
								  << a << ", " << b << ')'
								  << std::endl;
						std::cerr << __FILE__ << ':' << __LINE__
								  << ": information: "
								  << "f(a) = " << fa
								  << std::endl;
						std::cerr << __FILE__ << ':' << __LINE__
								  << ": information: "
								  << "b - a = " << b - a
								  << ", tolerance = " << tolerance
								  << std::endl;
					}
					return a;
				}
			}
		}
	}

	double
	solve_1d_aps(
		double (*function)(const double,
						   const std::vector<double_complex_t> &),
		const std::vector<double_complex_t> &coefficient,
		const double a, const double b, const double epsilon)
	{
		// Make a modifiable copy of the bounds
		double a_modified = a;
		double b_modified = b;
		// Constrain the epsilon to machine precision
		const double machine_epsilon =
			std::max(epsilon, machine_limit_t::spacing(a));
		const bool verbose = false;
		unsigned int niteration;
		unsigned int nevaluation;
		// Call solve_1d_aps_impl_rroot() to find the root
		const double root = solve_1d_aps_impl_rroot(
			function, coefficient, a_modified, b_modified,
			machine_epsilon, niteration, nevaluation, verbose);

		if(verbose) {
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": information: " << niteration
					  << " iteration(s), " << nevaluation
					  << " evaluations" << std::endl;

			const double function_value =
				function(root, coefficient);
			const double function_value_tolerance =
				4 * machine_limit_t::spacing(root);

			if(fabs(function_value) >
			   100 * function_value_tolerance) {
				std::cerr << __FILE__ << ':' << __LINE__
						  << ": warning: function value = "
						  << function_value << " > tolerance = "
						  << function_value_tolerance << std::endl;
			}
		}

		return root;
	}

}
