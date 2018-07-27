#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H
#include <jetbase/specfunc.h>

namespace jet {

#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 1572)
#endif // __INTEL_COMPILER

	// Special functions for kinematics

	double angular_range_reduce(const double x)
	{
		if(!std::isfinite(x))
			return x;

		I(F(x));

		static const double cody_waite_x_max = 1608.4954386379741381;
		// WARNING for PHENIX: DO NOT change two_pi_0 because it is
		// nearly 2 pi. The mantissa trunction IS INTENDED and needed
		// to allow multiplication against a small integer to be
		// exact. See W. Cody and W. Waite, Software Manual for the
		// Elementary Functions (Prentice-Hall, Englewood Cliffs, NJ,
		// 1980); W. J. Cody, IEEE Trans. Comput. C-22:6, 598--601
		// (1973).
		static const double two_pi_0 = 6.2831853071795649157;
		static const double two_pi_1 = 2.1561211432631314669e-14;
		static const double two_pi_2 = 1.1615423895917441336e-27;
		double ret;

		if(x >= -cody_waite_x_max && x <= cody_waite_x_max) {
			static const double inverse_two_pi =
				0.15915494309189534197;
			const double k = rint(x * inverse_two_pi);
			ret = ((x - (k * two_pi_0)) - k * two_pi_1) -
				k * two_pi_2;
		}
		else {
			long double sin_x;
			long double cos_x;

			sincosl(x, &sin_x, &cos_x);
			ret = (double)atan2l(sin_x, cos_x);
		}
		if(ret == -M_PI)
			ret = M_PI;

		if(!(ret > -M_PI && ret <= M_PI))
			fprintf(stderr, "ret = %.20f\n", ret);
		I(ret > -M_PI && ret <= M_PI);
		IG(FEQ(ret, x), x > -M_PI && x <= M_PI);
#ifndef NVERIFY
		I(FEQ(sinf(ret), sinf(x)));
		I(FEQ(cosf(ret), cosf(x)));
#endif // NVERIFY

		return ret;
	}

	float angular_range_reduce(const float x)
	{
		if(!std::isfinite(x))
			return x;

		I(F(x));

		static const float cody_waite_x_max = 1608.4954386F;
		// WARNING for PHENIX: DO NOT change two_pi_0 because it is
		// nearly 2 pi. The mantissa trunction IS INTENDED and needed
		// to allow multiplication against a small integer to be
		// exact. See W. Cody and W. Waite, Software Manual for the
		// Elementary Functions (Prentice-Hall, Englewood Cliffs, NJ,
		// 1980); W. J. Cody, IEEE Trans. Comput. C-22:6, 598--601
		// (1973).
		static const float two_pi_0 =  6.283203125F;
		static const float two_pi_1 = -1.781759784e-5F;
		static const float two_pi_2 = -2.22577267e-10F;
		float ret;

		if(x >= -cody_waite_x_max && x <= cody_waite_x_max) {
			static const float inverse_two_pi = 0.1591549431F;

			const float k = rintf(x * inverse_two_pi);
			ret = ((x - (k * two_pi_0)) - k * two_pi_1) -
				k * two_pi_2;
		}
		else {
			long double sin_x;
			long double cos_x;

			sincosl(x, &sin_x, &cos_x);
			ret = (float)atan2l(sin_x, cos_x);
		}
		if(ret == -(float)M_PI)
			ret = (float)M_PI;

		I(ret > -(float)M_PI && ret <= (float)M_PI);
		IG(FEQ(ret, x), x > -(float)M_PI && x <= (float)M_PI);
#ifndef NVERIFY
		I(FEQ(sin(ret), sin(x)));
		I(FEQ(cos(ret), cos(x)));
#endif // NVERIFY

		return ret;
	}

#ifdef _HALF_H_
	half angular_range_reduce(const half x)
	{
		// Arithmetics with IEEE 754r half precision bears no
		// advantage on ordinary CPU and is therefore pointless.
		return angular_range_reduce(float(x));
	}
#endif // _HALF_H_

#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif // __INTEL_COMPILER

	double hav(const double x)
	{
		const double sin_x_2 = sin(0.5 * x);
		const double retval = sin_x_2 * sin_x_2;

		return retval;
	}

	float havf(const float x)
	{
		const float sin_x_2 = sinf(0.5F * x);
		const float retval = sin_x_2 * sin_x_2;

		return retval;
	}

	double central_angle(const double l1, const double p1,
						 const double l2, const double p2)
	{
		const double hav_dl = hav(l2 - l1);
		const double hav_dp = hav(p2 - p1);
		const double cos_p1_cos_p2 = cos(p1) * cos(p2);
		const double retval = 2.0 *
			asin(sqrt(hav_dp + cos_p1_cos_p2 * hav_dl));

		return retval;
	}

	float central_angle(const float l1, const float p1,
						 const float l2, const float p2)
	{
		const float hav_dl = havf(l2 - l1);
		const float hav_dp = havf(p2 - p1);
		const float cos_p1_cos_p2 = cosf(p1) * cosf(p2);
		const float retval = 2.0F *
			asinf(sqrtf(hav_dp + cos_p1_cos_p2 * hav_dl));

		return retval;
	}

	double evaluate_chebyshev_series_nd_sparse(
		const std::vector<double> &x,
		const std::vector<double> &coefficient,
		const std::vector<std::vector<unsigned long> > &degree)
	{
		std::vector<std::vector<double> >
			t(x.size(), std::vector<double>(1, 1));

		for(unsigned long i = 0; i < x.size(); i++) {
			t[i].push_back(x[i]);
		}

		std::vector<double>::const_iterator coefficient_iterator =
			coefficient.begin();
		std::vector<std::vector<unsigned long> >::const_iterator
			degree_iterator = degree.begin();
		double sum = 0;

		for(; coefficient_iterator != coefficient.end() &&
				degree_iterator != degree.end();
			coefficient_iterator++,
				degree_iterator++) {
			double product = 1;

			I(degree_iterator->size() == x.size());

			for(unsigned long dimension = 0;
				dimension < degree_iterator->size(); dimension++) {
				const unsigned long next_degree =
					t[dimension].size();

				I(next_degree >= 2);

				const unsigned long n =
					(*degree_iterator)[dimension];

				// Iterative evaluation for T_2 and beyond using the
				// recurrent relation (e.g. Abramowith & Stegun,
				// 22.7.4):
				//
				// T_{n + 1}(x) = 2 x T_n(x) - T_{n - 1}(x), n >= 1
				if(n >= next_degree) {
					t[dimension].resize(n + 1);
					for(unsigned long i = next_degree; i <= n; i++) {
						t[dimension][i] =
							2 * x[dimension] * t[dimension][i - 1] -
							t[dimension][i - 2];
					}
				}
				product *= t[dimension][n];
			}
			sum += *coefficient_iterator * product;
		}

		return sum;
	}

	void evaluate_chebyshev_series_nd_sparse(
		double &f, double gradient[], double hessian[],
		const std::vector<double> &x,
		const std::vector<double> &coefficient,
		const std::vector<std::vector<unsigned long> > &degree)
	{
		std::vector<std::vector<double> >
			t(x.size(), std::vector<double>(1, 1));
		std::vector<std::vector<double> >
			u(x.size(), std::vector<double>(1, 1));
		std::vector<double> inverse_x_square_minus_1;

		inverse_x_square_minus_1.resize(x.size());
		for(unsigned long i = 0; i < x.size(); i++) {
			t[i].reserve(32);
			u[i].reserve(32);
			t[i].push_back(x[i]);
			u[i].push_back(2 * x[i]);
			inverse_x_square_minus_1[i] = 1 / (x[i] * x[i] - 1);
		}

		std::vector<double>::const_iterator
			coefficient_iterator = coefficient.begin();
		std::vector<std::vector<unsigned long> >::const_iterator
			degree_iterator = degree.begin();

		f = 0;
		for(unsigned long i = 0; i < x.size(); i++) {
			gradient[i] = 0;
		}
		for(unsigned long i = 0;
			i < ((x.size() * (x.size() + 1)) >> 1); i++) {
			hessian[i] = 0;
		}

		I(coefficient.size() == degree.size());

		for(; coefficient_iterator != coefficient.end() &&
				degree_iterator != degree.end();
			coefficient_iterator++,
				degree_iterator++) {
			const unsigned long dimension = degree_iterator->size();

			I(dimension == x.size());

			const unsigned long triangular_size =
				(dimension * (dimension + 1)) >> 1;
			double product_t = 1;
			std::vector<double> product_u(dimension, 1);
			std::vector<double> product_tu(triangular_size, 1);

			for(unsigned long i = 0; i < dimension; i++) {
				const unsigned long next_degree = t[i].size();

				I(next_degree >= 2);
				IG(u[i].size() == next_degree - 1, next_degree > 2);

				const unsigned long n = (*degree_iterator)[i];

				// Iterative evaluation for T_2 and beyond using the
				// recurrent relation (e.g. Abramowith & Stegun,
				// 22.7.4):
				//
				// T_{n + 1}(x) = 2 x T_n(x) - T_{n - 1}(x), n >= 1
				// U_{n + 1}(x) = 2 x U_n(x) - U_{n - 1}(x), n >= 1
				if(n >= next_degree) {
					t[i].resize(n + 1);
					for(unsigned long j = next_degree; j <= n; j++) {
						t[i][j] = 2 * x[i] *
							t[i][j - 1] - t[i][j - 2];
					}
					if(n >= 3) {
						u[i].resize(n);
						for(unsigned long j = next_degree - 1;
							j <= n - 1; j++) {
							u[i][j] = 2 * x[i] *
								u[i][j - 1] - u[i][j - 2];
						}
					}
				}
				const double t_n = t[i][n];
				const double n_u_n = n == 0 ? 0 : n * u[i][n - 1];
				const double n_n_t_minus_x_u_x_square_1 =
					n == 0 ? 0 :
					n * (n * t[i][n] - x[i] * u[i][n - 1]) *
					inverse_x_square_minus_1[i];

				// T
				product_t *= t_n;
				// U, T U diagonal
				for(unsigned long j = 0; j < dimension; j++) {
					product_u[j] *= j == i ? n_u_n : t_n;
					product_tu[j] *= j == i ?
						n_n_t_minus_x_u_x_square_1 : t_n;
				}
				// T U off diagonal
				for(unsigned long j = 1; j < dimension; j++) {
					for(unsigned long k = 0; k < j; k++) {
						const unsigned long index = dimension +
							((j * (j - 1)) >> 1) + k;

						I(index >= dimension);

						product_tu[index] *=
							(j == i || k == i) ? n_u_n : t_n;
					}
				}
			}
			f += *coefficient_iterator * product_t;
			for(unsigned long i = 0; i < dimension; i++)
				gradient[i] += *coefficient_iterator * product_u[i];
			for(unsigned long i = 0; i < triangular_size; i++)
				hessian[i] += *coefficient_iterator * product_tu[i];
		}
	}
}
