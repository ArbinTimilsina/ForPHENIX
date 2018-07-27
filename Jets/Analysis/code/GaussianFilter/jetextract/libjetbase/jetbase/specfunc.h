// -*- mode: c++; -*-

#ifndef JETBASE_SPECFUNC_H_
#define JETBASE_SPECFUNC_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <jetbase/dbc.h>
#include <jetbase/num.h>

// JET RECONSTRUCTION NUMERICS - SPECIAL FUNCTIONS

namespace jet {

	/**
	 * Returns the haversine of x
	 *
	 * @param[in] x argument
	 * @return the haversine of x, hav(x)
	 */
	double hav(const double x);
	/**
	 * Returns the haversine of x
	 *
	 * @param[in] x argument
	 * @return the haversine of x, hav(x)
	 */
	float havf(const float x);
	double central_angle(const double l1, const double p1,
						 const double l2, const double p2);
	float central_angle(const float l1, const float p1,
						const float l2, const float p2);

	// Kinematic functions

	/**
	 * Returns the angular range reduction of angle x
	 *
	 * @param[in] x angle
	 * @return the reduced angle y in (-pi, pi] such that y + 2 pi n =
	 * x
	 */
	double angular_range_reduce(const double x);
	/**
	 * Returns the angular range reduction of angle x
	 *
	 * @param[in] x angle
	 * @return the reduced angle y in (-pi, pi] such that y + 2 pi n =
	 * x
	 */
	float angular_range_reduce(const float x);

	double evaluate_chebyshev_series_nd_sparse(
		const std::vector<double> &x,
		const std::vector<double> &coefficient,
		const std::vector<std::vector<unsigned long> > &degree);
	void evaluate_chebyshev_series_nd_sparse(
		double &f, double gradient[], double hessian[],
		const std::vector<double> &x,
		const std::vector<double> &coefficient,
		const std::vector<std::vector<unsigned long> > &degree);
}

#endif // JETBASE_SPECFUNC_H_
