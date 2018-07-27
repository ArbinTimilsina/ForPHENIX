// -*- mode: c++; -*-

#ifndef XJETREC_SERIES_H_
#define XJETREC_SERIES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <jetbase/specfunc.h>

/////////////////////////////////////////////////////////////////////

namespace jet {

	class chebyshev_series_nd_sparse_t {
	private:
		std::vector<std::pair<double, double> > _range;
		std::vector<double> _coefficient;
		std::vector<std::vector<unsigned long> > _degree;
	public:
		inline chebyshev_series_nd_sparse_t(void)
		{
			I(empty());
		}
		inline chebyshev_series_nd_sparse_t(
			const std::vector<std::pair<double, double> > &range)
			: _range(range)
		{
		}
		inline std::vector<std::pair<double, double> > range() const
		{
			return _range;
		}
		inline void push_back(
			const double coefficient,
			const std::vector<unsigned long> &degree)
		{
			_coefficient.push_back(coefficient);
			_degree.push_back(degree);
		}
		inline bool empty(void) const
		{
			return _range.empty() || _coefficient.empty() ||
				_degree.empty();
		}
		inline void clear(void)
		{
			_coefficient.clear();
			_degree.clear();
		}
		inline double operator()(const std::vector<double> x) const
		{
			std::vector<double> x_scaled;

			x_scaled.resize(std::min(x.size(), _range.size()));
			if(x_scaled.size() == 0) {
				return NAN;
			}
			for(unsigned long i = 0; i < _range.size(); i++) {
				x_scaled[i] = 2 * (x[i] - _range[i].first) /
					(_range[i].second - _range[i].first) - 1;
			}

			return evaluate_chebyshev_series_nd_sparse(
				x_scaled, _coefficient, _degree);
		}
		inline void evaluate_order_2(
			double &f, double gradient[], double hessian[],
			const std::vector<double> x) const
		{
			std::vector<double> x_scaled;
			const unsigned long dimension = x.size();

			x_scaled.resize(dimension);
			for(unsigned long i = 0; i < dimension; i++)
				x_scaled[i] = 2 * (x[i] - _range[i].first) /
					(_range[i].second - _range[i].first) - 1;

			evaluate_chebyshev_series_nd_sparse(
				f, gradient, hessian,
				x_scaled, _coefficient, _degree);
			for(unsigned long i = 0; i < dimension; i++) {
				const double scale =
					2 / (_range[i].second - _range[i].first);

				gradient[i] *= scale;
				hessian[i] *= scale * scale;
			}
			for(unsigned long i = 1; i < dimension; i++)
				for(unsigned long j = 0; j < i; j++) {
					const unsigned long index =
						dimension + ((i * (i - 1)) >> 1) + j;

					hessian[index] *= 4 /
						((_range[i].second - _range[i].first) *
						 (_range[j].second - _range[j].first));
				}
		}
	};

}

#endif // XJETREC_SERIES_H_
