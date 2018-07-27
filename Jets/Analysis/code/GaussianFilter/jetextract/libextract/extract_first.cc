#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <extract/extract_first.h>
#include <jetstat/histogram.h>

/////////////////////////////////////////////////////////////////////

std::vector<jet::cluster_t> g_cluster_raw;

event_writer_first_t::event_writer_first_t(
	const bool legacy, const char *output_filename,
	const char *statistics_output_filename,
	const int collider_run_number,
	const jet::collisional_system_t collisional_system,
	const char *trigger, const unsigned int mode)
	: SubsysReco("xef_writer_first"),
	  _legacy(legacy), _nevent(0),
	  _collider_run_number(collider_run_number),
	  _experiment_run_number(-1),
	  _collisional_system(collisional_system),
	  _trigger(trigger), _output_filename(output_filename),
	  _output_file(NULL),
	  _mode(mode),
	  _statistics_output_filename(statistics_output_filename),
	  _statistics_output_file(NULL)
{
	_output_file = new jet::xef_file_t(
		_output_filename, _collider_run_number,
		_experiment_run_number, _collisional_system,
		_trigger, "Yue Shi Lai",
		"Columbia University and Nevis Laboratories",
		jet::xef_file_t::COMPRESSION_LZO1X);
}

event_writer_first_t::~event_writer_first_t(void)
{
}

int event_writer_first_t::Init(PHCompositeNode *top_node)
{
	return 0;
}

int event_writer_first_t::InitRun(PHCompositeNode *top_node)
{
	const RunHeader *phenix_run_header =
		findNode::getClass<RunHeader>(top_node, "RunHeader");

	if(phenix_run_header == NULL) {
		std::cerr << __FILE__ << ':' << __LINE__
				  << ": phenix_run_header == NULL" << std::endl;
		std::cerr << ">>> Begin: " << __FUNCTION__ << std::endl;
		return 1;
	}
	_experiment_run_number = phenix_run_header->get_RunNumber();
	std::cerr << __FILE__ << ':' << __LINE__
			  << ": _experiment_run_number = "
			  << _experiment_run_number << std::endl;
	std::cerr << __FILE__ << ':' << __LINE__
			  << ": _mode = " << _mode << std::endl;

	return 0;
}

int event_writer_first_t::process_event(PHCompositeNode *top_node)
{
	const PHGlobal *phenix_global =
		findNode::getClass<PHGlobal>(top_node, "PHGlobal");
	float bbc_t0 = NAN;

	if(phenix_global != NULL) {
		bbc_t0 = phenix_global->getBbcTimeZero();
	}

	emcClusterContainer *phenix_emc;

	if(_legacy) {
		PHTypedNodeIterator<emcClusterContainer> emc_iter(top_node);
		PHIODataNode<emcClusterContainer> *container =
			emc_iter.find("emcClusterContainer");

		phenix_emc = container != NULL ? container->getData() : NULL;
	}
	else {
		phenix_emc = findNode::getClass<emcClusterContainer>(
			top_node, "PhPhotonList");
	}
	if(phenix_emc == NULL && _nevent < 20) {
		std::cerr << __FILE__ << ':' << __LINE__ << ": "
			"information: phenix_emc == NULL" << std::endl;
	}
	else if(phenix_emc != NULL && phenix_emc->size() == 0 &&
			_nevent < 20) {
		std::cerr << __FILE__ << ':' << __LINE__ << ": "
			"information: phenix_emc->size() == 0" << std::endl;
	}

	std::vector<jet::cluster_t> cluster;

	if(phenix_emc != NULL) {
		for(int cluster_idx = 0;
			cluster_idx < (int)phenix_emc->size(); cluster_idx++) {
			const emcClusterContent *content =
				phenix_emc->getCluster(cluster_idx);

			// Momentum arithmetics
			const float energy = content->ecore();

			I(energy >= 0);

			const double x = content->x();
			const double y = content->y();
			const double perp = sqrt(x * x + y * y);
			const float azimuth = (float)atan2(y, x);
			const double z = content->z();
			const float polar_angle = (float)atan2(perp, z);

			I(polar_angle >= 0 && polar_angle <= M_PI);

			const float pseudorapidity =
				-log(tan(0.5F * polar_angle));

			I(azimuth >= -M_PI && azimuth <= 2 * M_PI);

			ID(jet::tower_id_t central_tower_id(
				content->arm(), content->sector(), content->iypos(),
				content->izpos()));
			std::vector<jet::tower_id_t> tower_id;
			std::vector<float> partial_energy_sum;

			for(int tower = 0; tower < content->multiplicity();
				tower++) {
				tower_id.push_back(content->towerid(tower));
				partial_energy_sum.push_back(content->partesum(tower));
			}

			I(tower_id.size() > 0);
			I(tower_id[0] == central_tower_id);

			const float electromagnetic_chi_square = content->chi2();
			const float tof_timing_mean =
				content->tof() <= -9999 ? NAN :
				content->tof() - bbc_t0;
			const float tof_timing_dispersion =
				content->tofdisp() <= -9999 ? NAN :
				content->tofdisp();
			const float tof_timing_min =
				content->tofmin() <= -9999 ? NAN :
				content->tofmin() - bbc_t0;
			const float tof_timing_max =
				content->tofmax() <= -9999 ? NAN :
				content->tofmax() - bbc_t0;
			const float incident_angle = polar_angle;
			const float corrected_dispersion_y =
				content->corrdispy();
			const float corrected_dispersion_z =
				content->corrdispz();

			I(polar_angle >= 0 && polar_angle <= M_PI);

			jet::cluster_t cluster1(
				energy, pseudorapidity,
				jet::angular_range_reduce(azimuth),
				cluster_idx, tower_id, partial_energy_sum,
				electromagnetic_chi_square,
				jet::tof_timing_t(
					tof_timing_mean, tof_timing_dispersion,
					tof_timing_min, tof_timing_max),
				incident_angle, corrected_dispersion_y,
				corrected_dispersion_z);

			cluster.push_back(cluster1);
		}
	}

	g_cluster_raw = cluster;
	// std::cerr << __FILE__ << ':' << __LINE__ << ": g_cluster_raw.size() = " << g_cluster_raw.size() << std::endl;
	_nevent++;

	return 0;
}

int event_writer_first_t::End(PHCompositeNode *top_node)
{
	return 0;
}
