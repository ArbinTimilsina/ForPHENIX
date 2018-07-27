#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H

#include <jetrec/background.h>
#include <iostream>

namespace jet {

#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 869)
#endif // __INTEL_COMPILER

	unsigned int factorized_background_model_t::
	vertex_discretize(const collision_geometry_t &geometry)
	{
		if(geometry.vertex_height() < -25.0F) {
			return 0;
		}
		if(!(geometry.vertex_height() < 25.0F)) {
			return nbin_vertex - 1;
		}
		else {
			return static_cast<unsigned int>(
				(geometry.vertex_height() + 25.0F) *
				(nbin_vertex / 50.0F));
		}
	}

	unsigned int factorized_background_model_t::
	vertex_centrality_discretize(
		const collision_geometry_t &geometry)
	{
		return vertex_discretize(geometry) * 101 +
			static_cast<unsigned int>(geometry.centrality());
	}

	double factorized_background_model_t::operator()(
		const collision_geometry_t &geometry,
		const double pseudorapidity, const double azimuth) const
	{
		if(empty() || !geometry.in_range() ||
		   pseudorapidity <=
		   _position_dependence[0].range()[0].first ||
		   pseudorapidity >=
		   _position_dependence[0].range()[0].second) {
			return 0;
		}

		std::vector<double> position;

		position.push_back(pseudorapidity);
		position.push_back(jet::angular_range_reduce(
			azimuth - 0.5 * M_PI));

		return vertex_centrality_dependence(geometry) *
			position_dependence(geometry)(position);
	}

	void factorized_background_model_t::evaluate_order_2(
		double &f, double gradient[], double hessian[],
		const collision_geometry_t &geometry,
		const double pseudorapidity, const double azimuth) const
	{
		if(empty() || !geometry.in_range() ||
		   !(pseudorapidity >=
			 _position_dependence[0].range()[0].first &&
			 pseudorapidity <=
			 _position_dependence[0].range()[0].second)) {
			f = 0;
			gradient[0] = 0;
			gradient[1] = 0;
			hessian[0] = 0;
			hessian[1] = 0;
			hessian[2] = 0;
			return;
		}

		std::vector<double> position;

		position.push_back(pseudorapidity);
		position.push_back(jet::angular_range_reduce(
			azimuth - 0.5 * M_PI));

		position_dependence(geometry).evaluate_order_2(
			f, gradient, hessian, position);

		const double scale = vertex_centrality_dependence(geometry);

		f *= scale;
		gradient[0] *= scale;
		gradient[1] *= scale;
		hessian[0] *= scale;
		hessian[1] *= scale;
		hessian[2] *= scale;
	}


#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif // __INTEL_COMPILER

#if 0
	void general_background_model_t::discretize_reaction_plane(
		unsigned int &index_low, unsigned int &index_high,
		double &weight_low, double &weight_high,
		const double reaction_plane)
	{
		const double reaction_plane_corrected =
			(0.5 * (jet::angular_range_reduce(
					2.0 * reaction_plane) + M_PI)) *
			(nbin_reaction_plane / M_PI) - 0.5;
		const double reaction_plane_corrected_floor =
			floor(reaction_plane_corrected);
		const double reaction_plane_corrected_ceil =
			ceil(reaction_plane_corrected);

		index_low =
			(static_cast<int>(reaction_plane_corrected_floor) +
			 nbin_reaction_plane) % nbin_reaction_plane;
		index_high =
			(static_cast<int>(reaction_plane_corrected_ceil) +
			 nbin_reaction_plane) % nbin_reaction_plane;
		weight_high = reaction_plane_corrected -
			reaction_plane_corrected_floor;
		weight_low = reaction_plane_corrected_ceil -
			reaction_plane_corrected;

		const double sum_weight = weight_low + weight_high;

		if(std::fpclassify(sum_weight) == FP_ZERO) {
			weight_low = 0.5;
			weight_high = 0.5;
		}
		else {
			weight_low /= sum_weight;
			weight_high /= sum_weight;
		}
	}

	double general_background_model_t::operator()(
		const double centrality, const double reaction_plane,
		const double pseudorapidity, const double azimuth) const
	{
		if(!(centrality >= 0 && centrality < nbin_centrality) ||
		   !std::isfinite(reaction_plane)) {
			return 0;
		}

		const unsigned int discrete_centrality =
			static_cast<unsigned int>(rint(centrality));
		unsigned int discrete_reaction_plane_low;
		unsigned int discrete_reaction_plane_high;
		double weight_low;
		double weight_high;

		discretize_reaction_plane(
			discrete_reaction_plane_low,
			discrete_reaction_plane_high,
			weight_low, weight_high, reaction_plane);

		const unsigned long index_low =
			discrete_centrality * nbin_reaction_plane +
			discrete_reaction_plane_low;
		const unsigned long index_high =
			discrete_centrality * nbin_reaction_plane +
			discrete_reaction_plane_high;

		if(!(index_low < 101 * 32 &&
			 index_high < 101 * 32)) {
			return 0;
		}

		std::vector<double> position;

		position.push_back(pseudorapidity);
		position.push_back(jet::angular_range_reduce(
			azimuth - 0.5 * M_PI));

		if(!(pseudorapidity >=
			 _position_dependence[index_low].range()[0].first &&
			 pseudorapidity <=
			 _position_dependence[index_low].range()[0].second)) {
			return 0;
		}

		return weight_low * _position_dependence[index_low](position) +
			weight_high * _position_dependence[index_high](position);
	}

	void general_background_model_t::evaluate_order_2(
		double &f, double gradient[], double hessian[],
		const double centrality, const double reaction_plane,
		const double pseudorapidity, const double azimuth) const
	{
		if(!(centrality >= 0 && centrality < nbin_centrality) ||
		   !std::isfinite(reaction_plane)) {
			f = 0;
			gradient[0] = 0;
			gradient[1] = 0;
			hessian[0] = 0;
			hessian[1] = 0;
			hessian[2] = 0;
			return;
		}

		const unsigned int discrete_centrality =
			static_cast<unsigned int>(rint(centrality));
		unsigned int discrete_reaction_plane_low;
		unsigned int discrete_reaction_plane_high;
		double weight_low;
		double weight_high;

		discretize_reaction_plane(
			discrete_reaction_plane_low,
			discrete_reaction_plane_high,
			weight_low, weight_high, reaction_plane);

		const unsigned long index_low =
			discrete_centrality * nbin_reaction_plane +
			discrete_reaction_plane_low;
		const unsigned long index_high =
			discrete_centrality * nbin_reaction_plane +
			discrete_reaction_plane_high;

		if(!(index_low < 101 * 32 &&
			 index_high < 101 * 32)) {
			f = 0;
			gradient[0] = 0;
			gradient[1] = 0;
			hessian[0] = 0;
			hessian[1] = 0;
			hessian[2] = 0;
			return;
		}
		if(!(pseudorapidity >=
			 _position_dependence[index_low].range()[0].first &&
			 pseudorapidity <=
			 _position_dependence[index_low].range()[0].second)) {
			f = 0;
			gradient[0] = 0;
			gradient[1] = 0;
			hessian[0] = 0;
			hessian[1] = 0;
			hessian[2] = 0;
			return;
		}

		std::vector<double> position;

		position.push_back(pseudorapidity);
		position.push_back(jet::angular_range_reduce(
			azimuth - 0.5 * M_PI));

		double f_low;
		double f_high;
		double gradient_low[2] __attribute__ ((aligned(16)));
		double gradient_high[2] __attribute__ ((aligned(16)));
		double hessian_low[4] __attribute__ ((aligned(16)));
		double hessian_high[4] __attribute__ ((aligned(16)));

		_position_dependence[index_low].evaluate_order_2(
			f_low, gradient_low, hessian_low, position);
		_position_dependence[index_high].evaluate_order_2(
			f_high, gradient_high, hessian_high, position);

		f = weight_low * f_low + weight_high * f_high;
		gradient[0] = weight_low * gradient_low[0] +
			weight_high * gradient_high[0];
		gradient[1] = weight_low * gradient_low[1] +
			weight_high * gradient_high[1];

		hessian[0] = weight_low * hessian_low[0] +
			weight_high * hessian_high[0];
		hessian[1] = weight_low * hessian_low[1] +
			weight_high * hessian_high[1];
		hessian[2] = weight_low * hessian_low[2] +
			weight_high * hessian_high[2];
	}
#endif

}

