#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H

#include <jetrec/fr.h>

namespace jet {

	float inner_gaussian_feature(
		const jet::jet_t &jet,
		const std::vector<jet::track_t> &track,
		jet::snowmass_vector_t &center,
		const float radius_inner)
	{
		static const float radius = 0.3F;

		center = jet.momentum();

		float sum_inner = 0;
		float count_inner = 0;
		std::vector<jet::track_t>::const_iterator track_max =
			track.end();
		float weighted_perp_square_max = -FLT_MAX;

		for(std::vector<jet::track_t>::const_iterator iterator =
				track.begin();
			iterator != track.end(); iterator++) {
			const float perp_05 =
				std::max(0.0F, iterator->momentum().perp());
			const float perp_square = perp_05 * perp_05;
			const float distance_square =
				iterator->momentum().radial_distance_square(center);
			const float distance = sqrtf(distance_square);
			const float weight =
				expf((-1.0F / (2.0F * radius * radius)) *
					 distance_square);
			const float weight_inner =
				expf((-1.0F / (2.0F * radius_inner * radius_inner)) *
					 distance_square);

			if(std::isfinite(perp_square)) {
				sum_inner += weight_inner * perp_square;
				count_inner += weight_inner;
				if(weight * perp_square > weighted_perp_square_max) {
					weighted_perp_square_max = weight * perp_square;
					track_max = iterator;
				}
			}
		}

		const float mean_inner = sum_inner;

		if(track_max == track.end()) {
			return mean_inner;
		}

		float sum_inner_off_center = 0;
		float count_inner_off_center = 0;

		center = track_max->momentum();
		for(std::vector<jet::track_t>::const_iterator iterator =
				track.begin();
			iterator != track.end(); iterator++) {
			const float perp_05 =
				std::max(0.0F, iterator->momentum().perp());
			const float perp_square = perp_05 * perp_05;
			const float distance_square =
				iterator->momentum().radial_distance_square(center);
			const float distance = sqrtf(distance_square);
			const float weight_inner =
				expf((-1.0F / (2.0F * radius_inner * radius_inner)) *
					 distance_square);

			if(std::isfinite(perp_square)) {
				sum_inner_off_center += weight_inner * perp_square;
				count_inner_off_center += weight_inner;
			}
		}

		const float mean_inner_off_center = sum_inner_off_center;

		if(mean_inner_off_center <= mean_inner) {
			center = jet.momentum();
		}

		return std::max(mean_inner, mean_inner_off_center);
	}

	float inner_gaussian_feature(
		const jet::jet_t &jet,
		const std::vector<jet::track_t> &track,
		jet::snowmass_vector_t &center,
		const float radius_inner,
		const float exponent,
		const bool adaption)
	{
		static const float radius = 0.3F;

		center = jet.momentum();

		float sum_inner = 0;
		float count_inner = 0;
		std::vector<jet::track_t>::const_iterator track_max =
			track.end();
		float weighted_perp_power_max = -FLT_MAX;

		for(std::vector<jet::track_t>::const_iterator iterator =
				track.begin();
			iterator != track.end(); iterator++) {
			const float perp =
				std::max(0.0F, iterator->momentum().perp());
			const float perp_power = powf(perp, exponent);
			const float distance_square =
				iterator->momentum().radial_distance_square(center);
			const float distance = sqrtf(distance_square);
			const float weight =
				expf((-1.0F / (2.0F * radius * radius)) *
					 distance_square);
			const float weight_inner =
				expf((-1.0F / (2.0F * radius_inner * radius_inner)) *
					 distance_square);

			if(std::isfinite(perp_power)) {
				sum_inner += weight_inner * perp_power;
				count_inner += weight_inner;
				if(weight * perp_power > weighted_perp_power_max) {
					weighted_perp_power_max = weight * perp_power;
					track_max = iterator;
				}
			}
		}

		const float mean_inner = sum_inner;

		if(track_max == track.end() || !adaption) {
			return mean_inner;
		}

		float sum_inner_off_center = 0;
		float count_inner_off_center = 0;

		center = track_max->momentum();
		for(std::vector<jet::track_t>::const_iterator iterator =
				track.begin();
			iterator != track.end(); iterator++) {
			const float perp =
				std::max(0.0F, iterator->momentum().perp());
			const float perp_power = powf(perp, exponent);
			const float distance_square =
				iterator->momentum().radial_distance_square(center);
			const float distance = sqrtf(distance_square);
			const float weight_inner =
				expf((-1.0F / (2.0F * radius_inner * radius_inner)) *
					 distance_square);

			if(std::isfinite(perp_power)) {
				sum_inner_off_center += weight_inner * perp_power;
				count_inner_off_center += weight_inner;
			}
		}

		const float mean_inner_off_center = sum_inner_off_center;

		if(mean_inner_off_center <= mean_inner) {
			center = jet.momentum();
		}

		return std::max(mean_inner, mean_inner_off_center);
	}

	float sum_jt_feature(const jet::jet_t &jet)
	{
		float sum = 0.0F;
		const std::vector<std::vector<jet::track_t>::const_iterator>
			constituent_outer = jet.constituent();

		for(std::vector<std::vector<
				jet::track_t>::const_iterator>::const_iterator
				iterator_outer = constituent_outer.begin();
			iterator_outer != constituent_outer.end();
			iterator_outer++) {
			const std::vector<std::vector<jet::track_t>::
				const_iterator> constituent_inner =
				(*iterator_outer)->constituent();

			for(std::vector<std::vector<
					jet::track_t>::const_iterator>::const_iterator
					iterator_inner = constituent_inner.begin();
				iterator_inner != constituent_inner.end();
				iterator_inner++) {
				sum += (*iterator_inner)->momentum().perp() *
					sinf((*iterator_inner)->momentum().
						 radial_distance(jet.momentum()));
			}
		}

		return sum;
	}

}
