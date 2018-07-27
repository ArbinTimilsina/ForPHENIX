#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H

#include <iostream>
#include <half.h>
#include <towerstat/towerstat.h>

/////////////////////////////////////////////////////////////////////

int tower_statistics_t::Init(PHCompositeNode *top_node)
{
	if(!_output_filename_spectrum.empty()) {
		_output_file_spectrum =
			fopen(_output_filename_spectrum.c_str(), "w");
	}

	I(_output_file_spectrum != NULL);

	if(!_output_filename_electron.empty()) {
		_output_file_electron =
			fopen(_output_filename_electron.c_str(), "w");
	}

	I(_output_file_electron != NULL);

	if(!_output_filename_tracking.empty()) {
		_output_file_tracking =
			fopen(_output_filename_tracking.c_str(), "w");

		I(_output_file_tracking != NULL);
	}

	if(!_output_filename_tof_histogram.empty()) {
		_output_file_tof_histogram =
			fopen(_output_filename_tof_histogram.c_str(), "w");

		I(_output_file_tof_histogram != NULL);
	}

	if(!_output_filename_tof_electron_histogram.empty()) {
		_output_file_tof_electron_histogram =
			fopen(_output_filename_tof_electron_histogram.c_str(), "w");

		I(_output_file_tof_electron_histogram != NULL);
	}

	if(!_output_filename_tof_hadron_histogram.empty()) {
		_output_file_tof_hadron_histogram =
			fopen(_output_filename_tof_hadron_histogram.c_str(), "w");

		I(_output_file_tof_hadron_histogram != NULL);
	}

        if(!_output_filename_bbc.empty()) {
                _output_file_bbc =
                        fopen(_output_filename_bbc.c_str(), "w");

                I(_output_file_bbc != NULL);
        }

	std::vector<jet::histogram_t<float, uint32_t>::fixed_bin_t>
		tof_histogram_bin;

	tof_histogram_bin.push_back(
		jet::histogram_t<float, uint32_t>::fixed_bin_t(
			"", 4, 0.1, 0.5));

	tof_histogram_bin.push_back(
		jet::histogram_t<float, uint32_t>::fixed_bin_t(
			"", 24768, 0, 24768));

	tof_histogram_bin.push_back(
		jet::histogram_t<float, uint32_t>::fixed_bin_t(
			"", 128, -1, 2));

	tof_histogram_bin.push_back(
		jet::histogram_t<float, uint32_t>::fixed_bin_t(
			"", 1024, -64, 64));

	_tof_histogram = new jet::histogram_t<float, uint32_t>(
		"", tof_histogram_bin);

	std::vector<jet::histogram_t<float, uint32_t>::fixed_bin_t>
		tof_pid_histogram_bin;

	tof_pid_histogram_bin.push_back(
		jet::histogram_t<float, uint32_t>::fixed_bin_t(
			"", 24768, 0, 24768));

	tof_pid_histogram_bin.push_back(
		jet::histogram_t<float, uint32_t>::fixed_bin_t(
			"", 128, -1, 2));

	tof_pid_histogram_bin.push_back(
		jet::histogram_t<float, uint32_t>::fixed_bin_t(
			"", 1024, -64, 64));

	_tof_electron_histogram = new jet::histogram_t<float, uint32_t>(
		"", tof_pid_histogram_bin);

	_tof_hadron_histogram = new jet::histogram_t<float, uint32_t>(
		"", tof_pid_histogram_bin);

	return 0;
}

int tower_statistics_t::InitRun(PHCompositeNode *top_node)
{
	RunHeader *phenix_run_header =
		findNode::getClass<RunHeader>(top_node, "RunHeader");

	if(phenix_run_header == NULL) {
		std::cerr << __FILE__ << ':' << __LINE__
				  << ": phenix_run_header == NULL" << std::endl;
		std::cerr << ">>> Begin: " << __FUNCTION__ << std::endl;
		return 1;
	}
	top_node->print();

	return 0;
}

int tower_statistics_t::process_event(PHCompositeNode *top_node)
{
	EventHeader *phenix_header =
		findNode::getClass<EventHeader>(top_node, "EventHeader");

	if(phenix_header == NULL) {
		std::cerr << __FILE__ << ':' << __LINE__
				  << ": phenix_header = NULL" << std::endl;
		std::cerr << ">>> End: " << __FUNCTION__ << std::endl;
		return 1;
	}

	PHGlobal *phenix_global =
		findNode::getClass<PHGlobal>(top_node, "PHGlobal");

	if(phenix_global == NULL) {
		std::cerr << __FILE__ << ':' << __LINE__
				  << ": phenix_global = NULL" << std::endl;
		std::cerr << ">>> End: " << __FUNCTION__ << std::endl;
		return 1;
	}

	const float centrality = phenix_global->getCentrality();

	if(_output_file_bbc != NULL) {
		uint32_t uint32[1];

		uint32[0] = phenix_header->get_EvtSequence();
		fwrite(uint32, sizeof(uint32_t), 1, _output_file_bbc);

		float real32[15];

		real32[0] = phenix_global->getCentrality();

		real32[1] = phenix_global->getBbcZVertex();
		real32[2] = phenix_global->getBbcZVertexError();
		real32[3] = phenix_global->getZdcZVertex();
		real32[4] = phenix_global->getZdcZVertexError();

		real32[5] = phenix_global->getBbcChargeN();
		real32[6] = phenix_global->getBbcChargeS();
		real32[7] = phenix_global->getBbcTimeN();
		real32[8] = phenix_global->getBbcTimeS();
		real32[9] = phenix_global->getBbcTimeZero();

		real32[10] = phenix_global->getZdcEnergyN();
		real32[11] = phenix_global->getZdcEnergyS();
		real32[12] = phenix_global->getZdcTimeN();
		real32[13] = phenix_global->getZdcTimeS();
		real32[14] = phenix_global->getZdcTimeZero();

		fwrite(real32, sizeof(float), 15, _output_file_bbc);
	}

	emcClusterContainer *phenix_emc =
		findNode::getClass<emcClusterContainer>(
			top_node, "PhPhotonList");
	if(phenix_emc == NULL) {
		phenix_emc = findNode::getClass<emcClusterContainer>(
			top_node, "EvPhPhotonList");
	}

	PHTypedNodeIterator<emcClusterContainer> *emc_iter = NULL;

	if(phenix_emc == NULL) {
		emc_iter =
			new PHTypedNodeIterator<emcClusterContainer>(top_node);
		phenix_emc = emc_iter->find("emcClusterContainer")->getData();
	}
	if(phenix_emc == NULL) {
		std::cerr << __FILE__ << ':' << __LINE__
				  << ": phenix_emc = NULL" << std::endl;
		std::cerr << ">>> End: " << __FUNCTION__ << std::endl;
		return 1;
	}

	const PHCentralTrack *phenix_central =
		findNode::getClass<PHCentralTrack>(
			top_node, "PHCentralTrack");

	if (phenix_central == NULL) {
		phenix_central =
                findNode::getClass<PHCentralTrack>(
                        top_node, "TrPhCglList");
	}

	if(phenix_central == NULL) {
		std::cerr << __FILE__ << ':' << __LINE__
				  << ": phenix_central = NULL" << std::endl;
		std::cerr << ">>> End: " << __FUNCTION__ << std::endl;
		return 1;
	}

	std::vector<std::pair<double, double> > tracking_hit;
	uint32_t dc_multiplicity_02 = 0;
	uint32_t dc_multiplicity_20 = 0;

	for(int track_idx = 0;
		track_idx < (int)phenix_central->get_npart(); track_idx++) {
		const int quality =
			phenix_central->get_quality(track_idx);
		const double polar_angle =
			phenix_central->get_the0(track_idx);
		const half pc3_sigma_displacement_azimuth =
			phenix_central->get_pc3sdphi(track_idx);
		const half pc3_sigma_displacement_height =
			phenix_central->get_pc3sdz(track_idx);
		const half emc_sigma_displacement_azimuth =
			phenix_central->get_emcsdphi(track_idx);
		const half emc_sigma_displacement_height =
			phenix_central->get_emcsdz(track_idx);
		const int rich_ring_nnormal_area_phototube =
			phenix_central->get_n0(track_idx);
		const int cluster_id =
			phenix_central->get_emcid(track_idx);
		const float momentum_norm =
			phenix_central->get_mom(track_idx);
		const float perp =
			(float)(momentum_norm * sin(polar_angle));

		if(quality > 7) {
			const double dc_crossing_azimuth =
				phenix_central->get_phi(track_idx);
			const double pc1_crossing_height =
				phenix_central->get_zed(track_idx);

			tracking_hit.push_back(std::pair<double, double>(
				dc_crossing_azimuth, pc1_crossing_height));

			if (std::isfinite(polar_angle)) {
				if (perp >= 0.2F) {
					dc_multiplicity_02++;
				}
				if (perp >= 2.0F) {
					dc_multiplicity_20++;
				}
			}
		}

		if(quality > 7 && std::isfinite(polar_angle) &&
		   polar_angle != -9999 &&
		   rich_ring_nnormal_area_phototube > 0 &&
		   cluster_id >= 0 &&
		   cluster_id < (int)phenix_emc->size()) {
			const half rich_ring_nnormal_area_photoelectron =
				phenix_central->get_npe0(track_idx);
			const int rich_ring_nwide_area_phototube =
				phenix_central->get_n1(track_idx);
			const half rich_ring_nwide_area_photoelectron =
				phenix_central->get_npe1(track_idx);
			const half rich_ring_unitful_chi_square =
				phenix_central->get_chi2(track_idx);
			const half rich_ring_displacement =
				phenix_central->get_disp(track_idx);
			const emcClusterContent *content =
				phenix_emc->getCluster(cluster_id);

			if(content != NULL) {
				const float energy = content->ecore();
				const half electromagnetic_chi_square =
					content->chi2();
				const half photon_probability =
					content->prob_photon();
				const short central_tower_id =
					content->towerid(0);

				fwrite(&central_tower_id, sizeof(short), 1,
					   _output_file_electron);
				fwrite(&momentum_norm, sizeof(float), 1,
					   _output_file_electron);

				unsigned short rich_info[6];

				rich_info[0] = rich_ring_nnormal_area_phototube;
				rich_info[1] = rich_ring_nnormal_area_photoelectron.bits();
				rich_info[2] = rich_ring_nwide_area_phototube;
				rich_info[3] = rich_ring_nwide_area_photoelectron.bits();
				rich_info[4] = rich_ring_unitful_chi_square.bits();
				rich_info[5] = rich_ring_displacement.bits();

				fwrite(rich_info, sizeof(unsigned short), 6,
					   _output_file_electron);

				fwrite(&energy, sizeof(float), 1,
					   _output_file_electron);

				unsigned short pc3_info[2];

				pc3_info[0] = pc3_sigma_displacement_height.bits();
				pc3_info[1] = pc3_sigma_displacement_azimuth.bits();

				fwrite(pc3_info, sizeof(unsigned short), 2,
					   _output_file_electron);

				unsigned short emc_info[4];

				emc_info[0] = emc_sigma_displacement_height.bits();
				emc_info[1] = emc_sigma_displacement_azimuth.bits();
				emc_info[2] = electromagnetic_chi_square.bits();
				emc_info[3] = photon_probability.bits();

				fwrite(emc_info, sizeof(unsigned short), 4,
					   _output_file_electron);

#if 0
				fprintf(_output_file_electron,
						"%.8e %d %.8e %d %.8e %.8e %.8e "
						"%.8e %.8e %.8e "
						"%.8e %.8e "
						"%d %.8e %.8e\n",
						momentum_norm,
						rich_ring_nnormal_area_phototube,
						rich_ring_nnormal_area_photoelectron,
						rich_ring_nwide_area_phototube,
						rich_ring_nwide_area_photoelectron,
						rich_ring_unitful_chi_square,
						rich_ring_displacement,
						//
						energy,
						electromagnetic_chi_square,
						photon_probability,
						//
						pc3_sigma_displacement_azimuth,
						pc3_sigma_displacement_height,
						//
						central_tower_id,
						energy / momentum_norm,
						energy * energy -
						momentum_norm * momentum_norm);
#endif
			}
		}
	}
	fwrite(&dc_multiplicity_02, sizeof(uint32_t), 1, _output_file_bbc);
	fwrite(&dc_multiplicity_20, sizeof(uint32_t), 1, _output_file_bbc);

	if(_output_file_tracking != NULL) {
		for(unsigned int i = 0; i < tracking_hit.size(); i++) {
			for(unsigned int j = 0; j < i; j++) {
				float d[2];

				d[0] = (tracking_hit[j].second -
						tracking_hit[i].second);
				d[1] = jet::angular_range_reduce(
							tracking_hit[j].first -
							tracking_hit[i].first);

				fwrite(d, sizeof(float), 2, _output_file_tracking);
#if 0
				fprintf(_output_file_tracking,
						"%.8e %.8e\n",
						jet::angular_range_reduce(
							tracking_hit[j].first -
							tracking_hit[i].first),
						tracking_hit[j].second -
						tracking_hit[i].second);
#endif
			}
		}
	}

	if(!(centrality >= 0 && centrality < 60)) {
		for(int cluster_idx = 0;
			cluster_idx < (int)phenix_emc->size(); cluster_idx++) {
			emcClusterContent *content =
				phenix_emc->getCluster(cluster_idx);

			// Momentum arithmetics
			const float energy = content->ecore();

			if(energy >= 0) {
				const int central_tower_id = content->towerid(0);

				if(energy < 32) {
					I(content->multiplicity() > 0);

					const int bin_energy = (int)floor(energy * 4);
					if(bin_energy < 128) {
						const int index =
							central_tower_id * _nbin_energy +
							bin_energy;

						I(index >= 0 &&
						  index < _ntower * _nbin_energy);

						_count[index]++;
					}
				}
				else {
					const int index =
						central_tower_id * _nbin_energy +
						(_nbin_energy - 1);

					I(index >= 0 && index < _ntower * _nbin_energy);

					_count[index]++;
				}
			}
		}
	}

	std::vector<bool> track_match(phenix_emc->size(), false);
	std::vector<bool> cluster_electron(phenix_emc->size(), false);
	std::vector<bool> cluster_hadron(phenix_emc->size(), false);

	for(int track_idx = 0;
		track_idx < (int)phenix_central->get_npart(); track_idx++) {
		const int quality =
			phenix_central->get_quality(track_idx);
		const double polar_angle =
			phenix_central->get_the0(track_idx);

		if(quality > 7 && std::isfinite(polar_angle) &&
		   polar_angle != -9999) {
			const int cluster_id =
				phenix_central->get_emcid(track_idx);

			if(cluster_id >= 0 && cluster_id <
			   static_cast<int>(phenix_emc->size())) {
				track_match[cluster_id] = true;

				const float track_i_pc3_dphi =
					phenix_central->get_pc3sdphi(track_idx);
				const float track_i_pc3_dz =
					phenix_central->get_pc3sdz(track_idx);
				const float track_i_pc3_displacement =
					sqrtf(track_i_pc3_dphi * track_i_pc3_dphi +
						  track_i_pc3_dz * track_i_pc3_dz);
				const bool track_i_is_electron =
					phenix_central->get_charge(track_idx) != 0 &&
					phenix_central->get_n0(track_idx) >=
					(track_i_pc3_displacement < 1.5F ? 2 : 1);

				if(track_i_is_electron) {
					cluster_electron[cluster_id] = true;
				}

				const float track_i_is_pion =
					phenix_central->get_isPi(track_idx) < 2.0F;
				const float track_i_is_kaon =
					phenix_central->get_isK(track_idx) < 2.0F;
				const float track_i_is_proton =
					phenix_central->get_isP(track_idx) < 2.0F;

				if(track_i_is_pion || track_i_is_kaon ||
				   track_i_is_proton) {
					cluster_hadron[cluster_id] = true;
				}
			}
		}
	}

	float bbc_t0 = NAN;

	if(phenix_global != NULL) {
		bbc_t0 = phenix_global->getBbcTimeZero();
	}

	for(int cluster_idx = 0;
		cluster_idx < (int)phenix_emc->size(); cluster_idx++) {
		const emcClusterContent *content =
			phenix_emc->getCluster(cluster_idx);
		const unsigned int central_tower_id = content->towerid(0);
		const float energy = content->ecore();
		const float tof_timing_mean = content->tof() <= -9999 ?
			NAN : content->tof() - bbc_t0;

		if(!track_match[cluster_idx]) {
			const float probability = content->prob_photon();
			std::vector<float> var(4U, NAN);

			var[0] = probability;
			var[1] = central_tower_id;
			var[2] = log10f(energy);
			var[3] = tof_timing_mean;

			reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
				_tof_histogram)->fill(var);
		}
		if(cluster_electron[cluster_idx]) {
			std::vector<float> var(3U, NAN);

			var[0] = central_tower_id;
			var[1] = log10f(energy);
			var[2] = tof_timing_mean;

			reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
				_tof_electron_histogram)->fill(var);
		}
		if(cluster_hadron[cluster_idx]) {
			std::vector<float> var(3U, NAN);

			var[0] = central_tower_id;
			var[1] = log10f(energy);
			var[2] = tof_timing_mean;

			reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
				_tof_hadron_histogram)->fill(var);
		}
	}

	if(emc_iter != NULL) {
		delete emc_iter;
	}

#if 0
	if(_nevent % (_nevent >= 1000 ? 1000 : _nevent >= 100 ? 100 :
				  _nevent >= 10 ? 10 : 1) == 0)
		std::cerr << "Run " << _experiment_run_number
				  << ": # of event : " << _nevent
				  << " event is processed" << std::endl;
#endif
	_nevent++;

	return 0;
}

int tower_statistics_t::End(PHCompositeNode *top_node)
{
	if(_output_file_spectrum != NULL) {
		fwrite(_count, sizeof(uint64_t), _ntower * _nbin_energy,
			   _output_file_spectrum);
		fclose(_output_file_spectrum);
		_output_file_spectrum = NULL;
	}
	if(_output_file_electron != NULL) {
		fclose(_output_file_electron);
		_output_file_electron = NULL;
	}
	if(_output_file_tracking != NULL) {
		fclose(_output_file_tracking);
		_output_file_tracking = NULL;
	}
	if(_output_file_tof_histogram != NULL) {
		jet::fwrite_histogram(
			*reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
				_tof_histogram),
			_output_file_tof_histogram);
		fclose(_output_file_tof_histogram);
		_output_file_tof_histogram = NULL;
		delete reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
			_tof_histogram);
		_tof_histogram = NULL;
	}
	if(_output_file_tof_electron_histogram != NULL) {
		jet::fwrite_histogram(
			*reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
				_tof_electron_histogram),
			_output_file_tof_electron_histogram);
		fclose(_output_file_tof_electron_histogram);
		_output_file_tof_electron_histogram = NULL;
		delete reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
			_tof_electron_histogram);
		_tof_electron_histogram = NULL;
	}
	if(_output_file_tof_hadron_histogram != NULL) {
		jet::fwrite_histogram(
			*reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
				_tof_hadron_histogram),
			_output_file_tof_hadron_histogram);
		fclose(_output_file_tof_hadron_histogram);
		_output_file_tof_hadron_histogram = NULL;
		delete reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
			_tof_hadron_histogram);
		_tof_hadron_histogram = NULL;
	}
        if(_output_file_bbc != NULL) {
                fclose(_output_file_bbc);
                _output_file_bbc = NULL;
        }


	return 0;
}
