// -*- mode: c++; -*-

#ifndef XJETREC_FR_H_
#define XJETREC_FR_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <cmath>
#include <vector>
#include <jetevent/particle.h>
#include <jetrec/jet.h>

namespace jet {

	extern float inner_gaussian_feature(
		const jet::jet_t &jet,
		const std::vector<jet::track_t> &track,
		jet::snowmass_vector_t &center,
		const float radius_inner = 0.1F);
	inline float inner_gaussian_feature(
		const jet::jet_t &jet,
		const std::vector<jet::track_t> &track,
		const float radius_inner = 0.1F)
	{
		jet::snowmass_vector_t center;

		return inner_gaussian_feature(
			jet, track, center, radius_inner);
	}
	extern float inner_gaussian_feature(
		const jet::jet_t &jet,
		const std::vector<jet::track_t> &track,
		jet::snowmass_vector_t &center,
		const float radius_inner,
		const float exponent,
		const bool adaption = true);
	inline float inner_gaussian_feature(
		const jet::jet_t &jet,
		const std::vector<jet::track_t> &track,
		const float radius_inner,
		const float exponent,
		const bool adaption = true)
	{
		jet::snowmass_vector_t center;

		return inner_gaussian_feature(
			jet, track, center, radius_inner,
			exponent, adaption);
	}
	extern float sum_jt_feature(const jet::jet_t &jet);

}

#endif // XJETREC_FR_H_
