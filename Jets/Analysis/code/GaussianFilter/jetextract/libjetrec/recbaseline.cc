#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600
#endif // _XOPEN_SOURCE
#include <cstdlib>
#include <functional>
#include <jetrec/rec.h>

namespace jet {

	std::vector<jet_t> reconstruction_leading_t::
	reconstruct(
		const std::vector<track_t>::const_iterator &track_begin,
		const std::vector<track_t>::const_iterator &track_end,
		const bool refine)
	{
		std::vector<track_t> track_sorted(track_begin, track_end);

		if(track_sorted.size() >= 2) {
			sort(track_sorted.begin(), track_sorted.end());
			track_sorted.resize(2);
		}

		std::vector<jet_t> ret;

		for(std::vector<track_t>::const_iterator
				iterator = track_sorted.begin();
			iterator != track_sorted.end(); iterator++) {
			if(iterator->momentum().perp() > FLT_EPSILON) {
				ret.push_back(iterator->momentum());
			}
		}

		return ret;
	} 

	std::vector<jet_t> reconstruction_lohsc_t::
	reconstruct(
		const std::vector<track_t>::const_iterator &track_begin,
		const std::vector<track_t>::const_iterator &track_end,
		const bool refine)
	{
		std::vector<jet_t> jet;
		std::vector<track_t> track_rest(track_begin, track_end);

		for(std::vector<track_t>::const_iterator
				iterator = track_begin;
			iterator != track_end; iterator++) {
			if(iterator->momentum().perp() > _seed_perp &&
			   fabsf(iterator->momentum().pseudorapidity()) <= 1.0F) {
				jet.push_back(iterator->momentum());
			}
		}
		std::sort(jet.begin(), jet.end());

		std::vector<jet_t> filtered_jet;

		for(std::vector<jet_t>::const_iterator iterator_outer =
				jet.begin();
			iterator_outer != jet.end(); iterator_outer++) {
			bool highest_in_cone = true;

			for(std::vector<jet_t>::const_iterator iterator_inner =
					jet.begin();
				iterator_inner != jet.end(); iterator_inner++) {
				if(iterator_inner->momentum().radial_distance(
						iterator_outer->momentum()) <
				   _inner_cone_radius &&
				   iterator_inner->momentum().perp() >
				   iterator_outer->momentum().perp()) {
					highest_in_cone = false;
					break;
				}
			}
			if(highest_in_cone) {
				filtered_jet.push_back(*iterator_outer);
			}
		}
		jet = filtered_jet;
		for(std::vector<jet_t>::iterator iterator_jet = jet.begin();
			iterator_jet != jet.end(); iterator_jet++) {
			float perp = 0.0F;
			std::vector<std::vector<track_t>::const_iterator> constituent;

			for(std::vector<track_t>::const_iterator iterator_track =
					track_begin;
				iterator_track != track_end; iterator_track++) {
				if(iterator_track->momentum().radial_distance(
					iterator_jet->momentum()) < _inner_cone_radius &&
				   iterator_track->momentum().perp() >=
				   _fragment_cut_perp) {
					constituent.push_back(iterator_track);
					perp += iterator_track->momentum().perp();
				}
			}
			iterator_jet->momentum().perp() = perp;
			iterator_jet->constituent() = constituent;
		}

		float perp_cone = 0;

		if(_background_subtract) {
			float perp_sum = 0.0F;
			unsigned long count_out_of_cone = 0;
			unsigned long count_all = 0;

			for(std::vector<track_t>::const_iterator iterator_track =
					track_begin;
				iterator_track != track_end; iterator_track++) {
				bool out_of_any_cone = true;

				for(std::vector<jet_t>::const_iterator iterator_jet =
						jet.begin();
					iterator_jet != jet.end(); iterator_jet++) {
					const bool in_cone = iterator_track->
						momentum().radial_distance(
							iterator_jet->momentum()) <
						_inner_cone_radius;

					if(in_cone) {
						out_of_any_cone = false;
						break;
					}
				}

				if(fabsf(iterator_track->momentum().
						 pseudorapidity()) < 1.0F) {
					if(out_of_any_cone) {
						count_out_of_cone++;
						if(iterator_track->momentum().perp() >=
						   _fragment_cut_perp) {
							perp_sum += iterator_track->momentum().perp();
						}
					}
					count_all++;
				}
			}

			const float perp_per_unit_area = perp_sum *
				(static_cast<float>(count_all) /
				 static_cast<float>(count_out_of_cone)) *
				static_cast<float>(1.0 / (4.0 * M_PI));

			perp_cone = perp_per_unit_area * M_PI *
				_inner_cone_radius * _inner_cone_radius;
		}
		if(_background_subtract) {
			fprintf(stderr, "%f\n", perp_cone);
			for(std::vector<jet_t>::iterator iterator = jet.begin();
				iterator != jet.end(); iterator++) {
				iterator->momentum().perp() -= perp_cone;
			}
		}

		return jet;
	}

}
