#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H
#include <stdint.h>
#include <jetevent/event.h>

namespace {

	static const unsigned int cut_statistics_nbin_cut = 64U;
	static const unsigned int cut_statistics_nbin_perp = 128U;
	static const float cut_statistics_log10_perp_min = -1.0F;
	static const float cut_statistics_log10_perp_max = 2.0F;
	static const float cut_statistics_log10_perp_scale =
		cut_statistics_nbin_perp /
		(cut_statistics_log10_perp_max -
		 cut_statistics_log10_perp_min);
	static const size_t _cut_statistics_size =
		cut_statistics_nbin_cut * (cut_statistics_nbin_perp + 2U);
	uint64_t _cut_statistics[_cut_statistics_size];

	template<typename p_t>
	inline void fill_cut_statistics(
		const std::vector<p_t> &particle,
		const unsigned int cut_id)
	{
		if(cut_id >= cut_statistics_nbin_cut) {
			return;
		}
		for(typename std::vector<p_t>::const_iterator iterator =
				particle.begin();
			iterator != particle.end(); iterator++) {
			const float log10_perp =
				log10f(iterator->momentum().perp());
			const unsigned int perp_bin =
				log10_perp < cut_statistics_log10_perp_min ?
				0U :
				!(log10_perp < cut_statistics_log10_perp_max) ?
				cut_statistics_nbin_perp + 1U :
				(log10_perp - cut_statistics_log10_perp_min) *
				cut_statistics_log10_perp_scale;
			const unsigned int index =
				cut_id * (cut_statistics_nbin_perp + 2U) + perp_bin;

			_cut_statistics[index]++;
		}
	}

	std::vector<jet::track_t> dc_matching_cut(
		bool &bad_event,
		const std::vector<jet::track_t> &track,
		const bool high_quality,
		const bool x1_x2_large_perp,
		const float max_magnitude = 25.0F)
	{
		std::vector<jet::track_t> retval;

		for(std::vector<jet::track_t>::const_iterator iterator =
				track.begin();
			iterator != track.end(); iterator++) {
			if(iterator->charge() != 0) {
				const float magnitude =
					iterator->momentum().cartesian_magnitude();
				const bool quality_ok = high_quality ?
					((iterator->quality() & 0x1f) == 0x1f) :
					(x1_x2_large_perp && magnitude >= 6.0F) ?
					((iterator->quality() & 0x38) != 0x0 &&
					 (iterator->quality() & 0x3) == 0x3) :
					((iterator->quality() & 0x38) != 0x0);

				if(quality_ok) {
					if(!(magnitude <= max_magnitude)) {
						bad_event = true;
#if 0
						// May be suppressed because this happens
						// rather frequently
						std::cerr << __FILE__ << ':' << __LINE__
								  << ": information: bad_event due "
							"to track with p = "
								  << magnitude << ", max_p = "
								  << max_magnitude << std::endl;
#endif
					}
					// Note that the PHENIX magic error number -9999
					// in IEEE 754-2008 half precision is -10000.
					else if(
						iterator->pc3_sigma_displacement().azimuth() >
						-9999.0F &&
						iterator->pc3_sigma_displacement().height() >
						-9999.0F) {
						const float pc3_displacement = iterator->
							pc3_sigma_displacement().radius();

						if(pc3_displacement <=
						   (magnitude >= 6.0F ? 2.5F : 3.0F)) {
							retval.push_back(*iterator);
						}
					}
				}
			}
			else {
				retval.push_back(*iterator);
			}
		}

		return retval;
	}

	std::vector<jet::track_t> dc_conversion_electron_pair_cut(
		const std::vector<jet::track_t> &track)
	{
		std::vector<bool> is_electron(track.size(), false);
		std::vector<bool> keep(track.size(), true);

		for(size_t i = 0; i < track.size(); i++) {
			const jet::track_t track_i = track[i];
			const float track_i_pc3_displacement =
				track_i.pc3_sigma_displacement().radius();
			const bool track_i_is_electron =
				track_i.charge() != 0 &&
				track_i.rich_ring().normal_area_count().
				nphototube() >=
				(track_i_pc3_displacement < 1.5F ? 2 : 1);

			if(track_i_is_electron) {
				for(size_t j = 0; j < i; j++) {
					const jet::track_t track_j = track[j];
					const float dpseudorapidity =
						fabsf(track_j.momentum().pseudorapidity() -
							  track_i.momentum().pseudorapidity());
					const float dazimuth =
						fabsf(jet::angular_range_reduce(
							track_j.momentum().azimuth() -
							track_i.momentum().azimuth()));
					const float track_j_pc3_displacement =
						track_j.pc3_sigma_displacement().radius();
					const bool track_j_is_electron =
						track_j.charge() != 0 &&
						track_j.rich_ring().normal_area_count().
						nphototube() >=
						(track_j_pc3_displacement < 1.5F ? 2 : 1);

					if(dpseudorapidity < 0.005F && dazimuth < 0.2F &&
					   track_j_is_electron) {
						keep[i] = false;
						keep[j] = false;
					}
				}
				if(track_i.momentum().perp() < 4.0F) {
					is_electron[i] = true;
				}
			}
		}

		std::vector<jet::track_t> retval;

		for(size_t i = 0; i < track.size(); i++) {
			if(keep[i]) {
				jet::track_t track_i = track[i];

				if(is_electron[i]) {
					track_i.pdg_code() = track_i.charge() > 0 ?
						jet::particle_t::PDG_CODE_POSITRON :
						jet::particle_t::PDG_CODE_ELECTRON;
				}
				retval.push_back(track_i);
			}
		}

		return retval;
	}

	std::vector<jet::track_t> dc_dead_region_cut(
		const std::vector<jet::track_t> &track)
	{
		static const size_t ndead_region = 37;
		static const float dead_region[ndead_region][4] = {
			{-0.5 * M_PI, -4, 1.5 * M_PI, 4},
			{-0.5 * M_PI, -90, 1.5 * M_PI, -78},
			{-0.5 * M_PI, 78, 1.5 * M_PI, 90},
			{-0.66, -100, -0.56, 100},
			{0.97, -100, 1.07, 100},
			{2.06, -100, 2.16, 100},
			{3.71, -100, 3.81, 100},
			{-0.39, -29, -0.34, -14},
			{-0.28, 0, -0.22, 10},
			{-0.21, -78, -0.19, 0},
			{-0.12, -78, -0.09, 78},
			{0.17, -78, 0.21, 78},
			{0.21, -78, 0.22, 0},
			{0.35, -27, 0.44, -15},
			{0.55, -27, 0.60, -15},
			{0.35, 62, 0.41, 73},
			{0.79, -63, 0.83, -51},
			{0.82, -51, 0.88, -21},
			{0.90, 0, 0.95, 10},
			{2.16, -78, 2.20, -71},
			{2.34, -72, 2.38, -45},
			{2.34, -24, 2.38, -12},
			{2.41, -55, 2.47, -42},
			{2.41, 0, 2.47, 18},
			{2.41, 64, 2.47, 73},
			{2.30, 0, 2.32, 78},
			{2.32, 0, 2.35, 18},
			{2.50, 18, 2.55, 29},
			{2.73, 0, 2.79, 9},
			{2.73, 52, 2.79, 78},
			{3.01, -78, 3.02, 0},
			{3.05, -72, 3.10, -58},
			{3.05, -27, 3.10, -15},
			{3.13, -78, 3.18, -70},
			{3.29, -78, 3.34, -70},
			{3.53, -10, 3.57, 0},
			{3.60, 35, 3.66, 45}
		};
		std::vector<bool> keep(track.size(), true);

		for(size_t i = 0; i < track.size(); i++) {
			const jet::track_t track_i = track[i];
			const float z = track_i.pc1_crossing_height();
			const float azimuth =
				track_i.dc_crossing().azimuth() <
				static_cast<float>(-0.5 * M_PI) ?
				track_i.dc_crossing().azimuth() +
				static_cast<float>(2.0 * M_PI) :
				track_i.dc_crossing().azimuth();

			for(size_t j = 0; j < ndead_region; j++) {
				if(azimuth >= dead_region[j][0] &&
				   azimuth < dead_region[j][2] &&
				   z >= dead_region[j][1] &&
				   z < dead_region[j][3]) {
					keep[i] = false;
					break;
				}
			}
		}

		std::vector<jet::track_t> retval;

		for(size_t i = 0; i < track.size(); i++) {
			if(keep[i]) {
				retval.push_back(track[i]);
			}
		}

		return retval;
	}

	std::vector<jet::track_t> dc_ghost_cut(
		const std::vector<jet::track_t> &track)
	{
		static const float inverse_dc_wire_distance =
			static_cast<float>(160.0 / M_PI);

		std::vector<bool> keep(track.size(), true);

		for(size_t i = 0; i < track.size(); i++) {
			const jet::track_t track_i = track[i];

			for(size_t j = 0; j < i; j++) {
				const jet::track_t track_j = track[j];
				const float dazimuth =
					fabsf(track_j.dc_crossing().azimuth() -
						  track_i.dc_crossing().azimuth()) *
					inverse_dc_wire_distance;
				const float dz =
					fabsf(track_j.pc1_crossing_height() -
						  track_i.pc1_crossing_height());

				// Acceleration for pairs outside the ghost cut region
				if(dz >= 1.0F || dazimuth >= 1.0F) {
					// Nothing
				}
				if((dz < 0.01F && dazimuth < 0.75F) ||
				   (dz < 0.02F &&
					dazimuth > 0.78F && dazimuth < 0.88F) ||
				   (dz < 0.075F && dazimuth < 0.35F) ||
				   (dz < 0.7F && dazimuth < 0.08F) ||
				   ((1.0F / (0.125F * 0.125F)) * dz * dz +
					(1.0F / (0.05F * 0.05F)) *
					(dazimuth - 0.5F) * (dazimuth - 0.5F) < 1.0F) || 
				   ((1.0F / (0.5F * 0.5F)) * dz * dz +
					(1.0F / (0.2F * 0.2F)) * dazimuth * dazimuth <
					1.0F)) {
					keep[i] = false;
				}
			}
		}

		std::vector<jet::track_t> retval;

		for(size_t i = 0; i < track.size(); i++) {
			if(keep[i]) {
				retval.push_back(track[i]);
			}
		}

		return retval;
	}

}

namespace jet {

	void reset_cut_statistics(void)
	{
		memset(_cut_statistics, 0, _cut_statistics_size);
	}

	const uint64_t *cut_statistics(void)
	{
		return _cut_statistics;
	}

	const size_t cut_statistics_size(void)
	{
		return _cut_statistics_size;
	}

	void event_t::apply_phenix_tracking_cut(
		bool &bad_event, const bool high_quality,
		const bool x1_x2_large_perp)
	{
		fill_cut_statistics(_track, 0U);
		_track = dc_matching_cut(bad_event, _track, high_quality,
								 x1_x2_large_perp, 25.0F);
		fill_cut_statistics(_track, 1U);
		_track = dc_conversion_electron_pair_cut(_track);
		fill_cut_statistics(_track, 2U);
		_track = dc_ghost_cut(_track);
		fill_cut_statistics(_track, 3U);
		_track = dc_dead_region_cut(_track);
	}

	std::vector<cluster_t> event_t::cluster_reconstruct(
		const bool tower_map[],
		const float cluster_tower_scale[]) const
	{
		static const float min_cluster_energy = 0.5F /* GeV */;
		static const float max_electromagnetic_chi_square = 3.0F;
		const size_t nraw_cluster = _cluster.size();
		std::vector<cluster_t> ret;

		// Convert remaining clusters to tracks
		for(size_t i = 0; i < nraw_cluster; i++) {
			if(_cluster[i].tower_id_match(tower_map)
#if 1
			   &&
			   _cluster[i].energy() > min_cluster_energy &&
			   _cluster[i].electromagnetic_chi_square() <
			   max_electromagnetic_chi_square
#endif
			   ) {
				snowmass_vector_t momentum =
					_cluster[i].energy_scale(
						cluster_tower_scale).momentum();

				ret.push_back(_cluster[i]);
				ret.back().momentum() = momentum;
			}
		}

		return ret;
	}

	std::vector<cluster_t> event_t::cluster_reconstruct_no_recal(void) const
	{
		static const float min_cluster_energy = 0.5F /* GeV */;
		static const float max_electromagnetic_chi_square = 3.0F;
		const size_t nraw_cluster = _cluster.size();
		std::vector<cluster_t> ret;

		// Convert remaining clusters to tracks
		for(size_t i = 0; i < nraw_cluster; i++) {
			snowmass_vector_t momentum =
				_cluster[i].momentum();

			ret.push_back(_cluster[i]);
			ret.back().momentum() = momentum;
		}

		return ret;
	}

	event_t event_t::phenix_reconstruct(
		bool &bad_event, const bool tower_map[],
		const float cluster_tower_scale[],
		const float run_timing_mean[],
		const float tower_timing_mean[],
		const float tower_timing_standard_deviation[],
		const bool high_quality, const bool x1_x2_large_perp,
		const bool apply_timing_cut,
		const bool reconstruct_neutral_pion) const
	{
		static const float min_cluster_energy = 0.5F /* GeV */;
		static const float max_electromagnetic_chi_square = 3.0F;
		static const std::pair<float, float> pi_0_mass_range =
			std::pair<float, float>(
				1.3533535e-01F - 3.0F * 9.8774734e-03,
				1.3533535e-01F + 3.0F * 9.8774734e-03) /* GeV/c^2 */;
		const size_t nraw_cluster = _cluster.size();
		event_t event = *this;

		bad_event = false;
		if(event.centrality() <= -9999) {
			event.centrality() = NAN;
		}
		if(!(fabsf(event.vertex_height().beam_beam_counter_mean())
			 < 25.0F)) {
			bad_event = true;
#if 0
			// May be suppressed because this happens rather
			// frequently
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": information: event_number = "
					  << event._event_number
					  << ": bad_event due to vertex_height = "
					  << event.vertex_height()
					  << std::endl;
#endif
		}

		event.apply_phenix_tracking_cut(
			bad_event, high_quality, x1_x2_large_perp);

		const size_t ntrack_1 = event._track.size();

		for(size_t j = 0; j < ntrack_1; j++) {
			for(size_t k = 0; k < track_t::nuser_data; k++) {
				event._track[j].user_data()[k] = LLONG_MIN;
			}
			if(event._track[j].charge() == 0) {
				event._track[j].charge() = INT_MIN;
			}
		}

		// PHENIX tracking/calorimetry matching

		std::vector<bool> track_match_cluster(ntrack_1, false);
		std::vector<bool> track_multiple_match(ntrack_1, false);
		// Helper variables for the first pass
		std::vector<bool> cluster_match_track(nraw_cluster, false);
		std::vector<bool> cluster_multiple_match(nraw_cluster, false);

		// Track to cluster forward matching, pass 1
		for(size_t i = 0; i < ntrack_1; i++) {
			int cluster_id = event._track[i].cluster_id();

			if(cluster_id >= 0 && cluster_id < (int)nraw_cluster) {
				if(cluster_match_track[cluster_id]) {
					cluster_multiple_match[cluster_id] = true;
				}
				cluster_match_track[cluster_id] = true;
				event._track[i].energy() =
					_cluster[cluster_id].energy();
				event._track[i].user_data()[0] =
					_cluster[cluster_id].central_tower_id();
				track_match_cluster[i] = true;
			}
		}
		// Track to cluster forward matching, pass 2
		for(size_t i = 0; i < ntrack_1; i++) {
			int cluster_id = event._track[i].cluster_id();

			if(cluster_id >= 0 && cluster_id < (int)nraw_cluster &&
			   cluster_multiple_match[cluster_id]) {
				track_multiple_match[i] = true;
			}
		}
		// Convert remaining clusters to tracks
		for(size_t i = 0; i < nraw_cluster; i++) {
			const float tof_timing_offset =
				(run_timing_mean == NULL ||
				 tower_timing_mean == NULL) ? NAN :
				run_timing_mean[
					_cluster[i].central_tower_id().sector()] +
				tower_timing_mean[
					static_cast<int>(_cluster[i].central_tower_id())];
			const float tof_timing_standard_deviation =
				tower_timing_standard_deviation == NULL ? NAN :
				tower_timing_standard_deviation[
					static_cast<int>(_cluster[i].central_tower_id())];
			const float tof_timing_limit_0 =
				tof_timing_offset +
				-3.0F * tof_timing_standard_deviation;
			const float tof_timing_limit_1 =
				tof_timing_offset +
				3.0F * tof_timing_standard_deviation +
				0.9F / _cluster[i].energy();

			if(!cluster_match_track[i] &&
			   _cluster[i].tower_id_match(tower_map) &&
			   _cluster[i].energy() > min_cluster_energy &&
			   _cluster[i].electromagnetic_chi_square() <
			   max_electromagnetic_chi_square &&
			   (!apply_timing_cut ||
				(_cluster[i].tof_timing().mean() >=
				 tof_timing_limit_0 &&
				 _cluster[i].tof_timing().mean() <
				 tof_timing_limit_1))) {
				snowmass_vector_t momentum =
					_cluster[i].energy_scale(
						cluster_tower_scale).momentum();

				momentum.set_lightlike_perp();

				track_t track = track_t(momentum);

				track.charge() = 0;
				track.pdg_code() = particle_t::PDG_CODE_PHOTON;
				track.status_code() =
					_cluster[i].status_code();
				track.user_data()[0] =
					_cluster[i].central_tower_id();

				float cluster_data[4];

				cluster_data[0] = _cluster[i].electromagnetic_chi_square();
				cluster_data[1] = _cluster[i].corrected_dispersion_y();
				cluster_data[2] = _cluster[i].corrected_dispersion_z();
				cluster_data[3] = _cluster[i].tof_timing().mean();

				track.user_data()[1] =
					*reinterpret_cast<uint32_t *>(cluster_data);
				track.user_data()[2] =
					*reinterpret_cast<uint32_t *>(cluster_data + 1);
				track.user_data()[3] =
					*reinterpret_cast<uint32_t *>(cluster_data + 2);
				track.user_data()[4] =
					*reinterpret_cast<uint32_t *>(cluster_data + 3);
				event._track.push_back(track);
			}
		}

		// Add all normal tracks

		const size_t ntrack_2 = event._track.size();
		std::vector<track_t> new_track;

		new_track.reserve(ntrack_2);
		for(size_t i = 0; i < ntrack_1; i++) {
			if(std::isfinite(
				event._track[i].momentum().pseudorapidity()) &&
			   std::isfinite(
				event._track[i].momentum().azimuth())) {
				bool valid_track;

				if(track_multiple_match[i]) {
					// Ambiguous track-cluster association
					event._track[i].energy() = NAN;
					valid_track = true;
				}
				else if(track_match_cluster[i]) {
					// E/p rejection of PHENIX conversion electrons
					if(std::isfinite(event._track[i].energy()) &&
					   event._track[i].rich_ring().
					   normal_area_count().nphototube() > 1 &&
					   event._track[i].momentum().energy() <
					   0.125F * event._track[i].momentum().
					   cartesian_magnitude()) {
						event._track[i].momentum().
							set_lightlike_perp();
					}
					valid_track = true;
				}
				else {
					// No EMC response
					event._track[i].momentum().set_lightlike_time();

					const float pc3_displacement =
						event._track[i].pc3_sigma_displacement().
						radius();
					const bool is_electron =
						event._track[i].rich_ring().
						normal_area_count().nphototube() >=
						(pc3_displacement < 1.5F ? 2 : 1);
					valid_track = !is_electron;
					if(is_electron) {
						bad_event = true;
#if 0
						// May be suppressed because this happens
						// rather frequently
						std::cerr << __FILE__ << ':' << __LINE__
								  << ": information: event_number = "
								  << event._event_number
								  << ": bad_event due to electron "
							"without EMC response"
								  << std::endl;
#endif
					}
				}
				if(valid_track) {
					new_track.push_back(event._track[i]);
				}
			}
		}
		fill_cut_statistics(new_track, 4U);
		for(size_t i = ntrack_1; i < ntrack_2; i++) {
			new_track.push_back(event._track[i]);
		}

		// Neutral pion reconstruction

		const size_t ntrack_3 = new_track.size();
		std::vector<bool> pi_0_match(ntrack_3, false);

		if(reconstruct_neutral_pion) {
			for(size_t i = 0; i < ntrack_3; i++) {
				for(size_t j = 0; j < i; j++) {
					if(new_track[i].pdg_code() ==
					   particle_t::PDG_CODE_PHOTON &&
					   new_track[j].pdg_code() ==
					   particle_t::PDG_CODE_PHOTON) {
						const lorentz_vector_t momentum_i =
							new_track[i].momentum();
						const lorentz_vector_t momentum_j =
							new_track[j].momentum();
						const lorentz_vector_t momentum_sum =
							momentum_i + momentum_j;
						const float mass = momentum_sum.mass();

						if(mass > pi_0_mass_range.first &&
						   mass < pi_0_mass_range.second &&
						   !pi_0_match[i] && !pi_0_match[j]) {
							pi_0_match[i] = true;
							pi_0_match[j] = true;

							track_t track = track_t(
								snowmass_vector_t(momentum_sum));

							track.charge() = 0;
							track.pdg_code() =
								particle_t::PDG_CODE_PION_0;
							track.status_code() =
								new_track[i].status_code() == 1 ?
								new_track[j].status_code() :
								new_track[i].status_code();
							new_track.push_back(track);
						}
					}
				}
			}
		}
		event._track.clear();
		event._cluster.clear();

		const size_t ntrack_4 = new_track.size();

		event._track.reserve(ntrack_4);
		for(size_t i = 0; i < ntrack_3; i++) {
			if(!pi_0_match[i]) {
				event._track.push_back(new_track[i]);
			}
		}
		for(size_t i = ntrack_3; i < ntrack_4; i++) {
			event._track.push_back(new_track[i]);
		}

		return event;
	}

}
