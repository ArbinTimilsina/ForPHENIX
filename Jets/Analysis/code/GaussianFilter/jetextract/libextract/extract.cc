#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <extract/extract.h>
#include <jetstat/histogram.h>

/////////////////////////////////////////////////////////////////////

extern std::vector<jet::cluster_t> g_cluster_raw;

namespace {

	float inner_gaussian_feature(
		const jet::jet_t &jet, const std::vector<jet::track_t> &track,
		const float radius_inner = 0.1F)
	{
		static const float radius = 0.3F;
		jet::snowmass_vector_t center = jet.momentum();

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
			const float weight_inner =
				expf((-1.0F / (2.0F * radius_inner * radius_inner)) *
					 distance_square);

			if(std::isfinite(perp_square)) {
				sum_inner_off_center += weight_inner * perp_square;
				count_inner_off_center += weight_inner;
			}
		}

		const float mean_inner_off_center = sum_inner_off_center;

		return std::max(mean_inner, mean_inner_off_center);
	}

	std::vector<jet::jet_t> gaussian_fake_rejection(
		const std::vector<jet::jet_t> &jet,
		const jet::event_t &event, const float threshold)
	{
		std::vector<jet::jet_t> retval;

		for(std::vector<jet::jet_t>::const_iterator iterator =
				jet.begin();
			iterator != jet.end(); iterator++) {
			if(inner_gaussian_feature(*iterator, event.track()) >=
			   threshold) {
				retval.push_back(*iterator);
			}
		}

		return retval;
	}

#if 1
	void set_run_5_cu_cu_background_generic(
		jet::factorized_background_model_t &model)
	{
#include <extract/run_5_cu_cu_background_centrality.h>

		for(unsigned int i = 0U; i < 101U * 10U; i++) {
			model.vertex_centrality_dependence()[i] =
				vertex_centrality_data[i];
		}

		std::vector<std::pair<double, double> > position_range;

		position_range.push_back(std::pair<double, double>(
			-0.525, 0.525));
		position_range.push_back(std::pair<double, double>(
			-M_PI, M_PI));

		for(unsigned int i = 0; i < 10; i++) {
			model.position_dependence()[i] = 
				jet::chebyshev_series_nd_sparse_t(position_range);
		}
	}

	void set_run_5_cu_cu_03_background(
		jet::factorized_background_model_t &model)
	{
		set_run_5_cu_cu_background_generic(model);

		std::vector<unsigned long> degree(2, 0);

#include <extract/run_5_cu_cu_background_position.h>

		for(unsigned int i = 0; i < 10; i++) {
			for(unsigned int j = 0; j < 151; j++) {
				const unsigned int index = i * 151 + j;

				degree[0] = vertex_position_data[index]._degree_0;
				degree[1] = vertex_position_data[index]._degree_1;
				model.position_dependence()[i].push_back(
					vertex_position_data[index]._x, degree);
			}
		}
	}

	void set_run_5_cu_cu_05_background(
		jet::factorized_background_model_t &model)
	{
	}
#endif

	void read_reattempt(void *array, const size_t size,
						const size_t nmemb, const char *filename)
	{
		int read_size = -1;
		FILE *fp = fopen(filename, "r");

		if(fp != NULL) {
			read_size = fread(array, size, nmemb, fp);
			fclose(fp);
			fp = NULL;
		}
		std::cerr << __FILE__ << ':' << __LINE__
				  << ": information: read_size = " << read_size
				  << ", nmemb = " << nmemb << std::endl;

		unsigned int reattempt = 0;

		while(read_size != static_cast<int>(nmemb)) {
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": information: reattempt in 5 sec"
					  << std::endl;
			sleep(5);
			fp = fopen(filename, "r");
			if(fp != NULL) {
				read_size = fread(array, size, nmemb, fp);
				fclose(fp);
				fp = NULL;
			}
			reattempt++;
			if(reattempt >= 64) {
				std::cerr << __FILE__ << ':' << __LINE__
						  << ": warning: reattempts exhausted"
						  << std::endl;
				exit(EXIT_FAILURE);
			}
		}
	}

	void set_run_7_au_au_03_background(
		jet::general_background_model_t &model)
	{
		static const unsigned long __polynomial_power[151 * 2] = {
			0,  0,
			0,  1,
			0,  2,
			0,  3,
			0,  4,
			0,  5,
			0,  6,
			0,  7,
			0,  8,
			0,  9,
			0, 10,
			0, 11,
			0, 12,
			0, 13,
			0, 14,
			0, 15,
			1,  0,
			1,  1,
			1,  2,
			1,  3,
			1,  4,
			1,  5,
			1,  6,
			1,  7,
			1,  8,
			1,  9,
			1, 10,
			1, 11,
			1, 12,
			1, 13,
			1, 14,
			1, 15,
			2,  0,
			2,  1,
			2,  2,
			2,  3,
			2,  4,
			2,  5,
			2,  6,
			2,  7,
			2,  8,
			2,  9,
			2, 10,
			2, 11,
			2, 12,
			2, 13,
			2, 14,
			3,  0,
			3,  1,
			3,  2,
			3,  3,
			3,  4,
			3,  5,
			3,  6,
			3,  7,
			3,  8,
			3,  9,
			3, 10,
			3, 11,
			3, 12,
			3, 13,
			4,  0,
			4,  1,
			4,  2,
			4,  3,
			4,  4,
			4,  5,
			4,  6,
			4,  7,
			4,  8,
			4,  9,
			4, 10,
			4, 11,
			4, 12,
			5,  0,
			5,  1,
			5,  2,
			5,  3,
			5,  4,
			5,  5,
			5,  6,
			5,  7,
			5,  8,
			5,  9,
			5, 10,
			5, 11,
			6,  0,
			6,  1,
			6,  2,
			6,  3,
			6,  4,
			6,  5,
			6,  6,
			6,  7,
			6,  8,
			6,  9,
			6, 10,
			7,  0,
			7,  1,
			7,  2,
			7,  3,
			7,  4,
			7,  5,
			7,  6,
			7,  7,
			7,  8,
			7,  9,
			8,  0,
			8,  1,
			8,  2,
			8,  3,
			8,  4,
			8,  5,
			8,  6,
			8,  7,
			8,  8,
			9,  0,
			9,  1,
			9,  2,
			9,  3,
			9,  4,
			9,  5,
			9,  6,
			9,  7,
			10,  0,
			10,  1,
			10,  2,
			10,  3,
			10,  4,
			10,  5,
			10,  6,
			11,  0,
			11,  1,
			11,  2,
			11,  3,
			11,  4,
			11,  5,
			12,  0,
			12,  1,
			12,  2,
			12,  3,
			12,  4,
			13,  0,
			13,  1,
			13,  2,
			13,  3,
			14,  0,
			14,  1,
			14,  2,
			15,  0,
			15,  1
		};
		static const unsigned int ncoefficient = 94 * 32 * 151;
		static const char *filename =
			"/phenix/hp/data10/ylai/train/data/"
			"run_7_au_au_background_coefficient.bin";
		double polynomial_coefficient[ncoefficient];

		read_reattempt(polynomial_coefficient, sizeof(double),
					   ncoefficient, filename);

		std::vector<std::pair<double, double> > position_range;

		position_range.push_back(std::pair<double, double>(
			-0.525, 0.525));
		position_range.push_back(std::pair<double, double>(
			-M_PI, M_PI));

		jet::chebyshev_series_nd_sparse_t position_dependence[101 * 32];

		for(int i = 0; i < 101 * 32; i++) {
			position_dependence[i] =
				jet::chebyshev_series_nd_sparse_t(position_range);

			std::vector<unsigned long> degree(2, 0);

			if(i < 32 || i >= 94 * 32) {
				degree[0] =  0;
				degree[1] =  0;
				position_dependence[i].push_back(0, degree);
			}
			else {
				for(int j = 0; j < 151; j++) {
					degree[0] = __polynomial_power[2 * j];
					degree[1] = __polynomial_power[2 * j + 1];
					position_dependence[i].push_back(
						polynomial_coefficient[(i - 32) * 151 + j],
						degree);
				}
			}
		}
		model = position_dependence;
	}
}



event_writer_t::event_writer_t(
	const bool legacy, const char *output_filename,
	const char *statistics_output_filename,
	const int collider_run_number,
	const jet::collisional_system_t collisional_system,
	const char *trigger, const unsigned int mode)
	: SubsysReco("xef_writer"),
	  _legacy(legacy), _nevent(0),
	  _collider_run_number(collider_run_number),
	  _experiment_run_number(-1),
	  _collisional_system(collisional_system),
	  _trigger(trigger), _output_filename(output_filename),
	  _output_file(NULL),
	  _mode(mode),
	  _filter_03(NULL), _filter_05(NULL),
	  _kt_03(NULL), _kt_05(NULL),
	  _antikt_03(NULL), _antikt_05(NULL),
	  _siscone_03(NULL), _siscone_05(NULL),
	  _statistics_output_filename(statistics_output_filename),
	  _statistics_output_file(NULL)
{
	_output_file = new jet::xef_file_t(
		_output_filename, _collider_run_number,
		_experiment_run_number, _collisional_system,
		_trigger, "Yue Shi Lai",
		"Columbia University and Nevis Laboratories",
		jet::xef_file_t::COMPRESSION_LZO1X);

	_background_model_cu_cu_03_perp =
			new jet::factorized_background_model_t();
	_background_model_cu_cu_05_perp =
			new jet::factorized_background_model_t();
	_background_model_au_au_03_perp =
			new jet::general_background_model_t();
	_background_model_03_time =
			new jet::factorized_background_model_t();
	_background_model_03_z =
			new jet::factorized_background_model_t();
	_background_model_05_time =
			new jet::factorized_background_model_t();
	_background_model_05_z =
			new jet::factorized_background_model_t();

	static const std::pair<float, float>
		phenix_acceptance_pseudorapidity_range =
		std::pair<float, float>(-0.525F, 0.525F);
	static const int phenix_acceptance_npixel_pseudorapidity = 42;

	switch(_mode) {
	case JET_RECONSTRUCT_RUN_5_CU_CU:
		set_run_5_cu_cu_03_background(*_background_model_cu_cu_03_perp);
		set_run_5_cu_cu_05_background(*_background_model_cu_cu_05_perp);
		break;
	case JET_RECONSTRUCT_RUN_7_AU_AU:
		set_run_7_au_au_03_background(*_background_model_au_au_03_perp);
		break;
	}

	if(_mode == JET_RECONSTRUCT_RUN_7_AU_AU) {
		static const char *filename =
			"/phenix/hp/data10/ylai/train/data/"
			"run_7_au_au_background_cache.bin";
		static const unsigned int npixel = 101 * 32 * 42 * 256;
		float *background_cache = new float[npixel];

		read_reattempt(background_cache, sizeof(float), npixel,
					   filename);
		_filter_03 = new jet::reconstruction_filtering_iir_t(
			0.3F, phenix_acceptance_pseudorapidity_range,
			phenix_acceptance_npixel_pseudorapidity,
			*_background_model_au_au_03_perp,
			*_background_model_03_time,
			*_background_model_03_z, background_cache);

		delete [] background_cache;
	}
	else {
		_filter_03 = new jet::reconstruction_filtering_iir_t(
			0.3F, phenix_acceptance_pseudorapidity_range,
			phenix_acceptance_npixel_pseudorapidity,
			*_background_model_cu_cu_03_perp,
			*_background_model_03_time,
			*_background_model_03_z);
	}
	_filter_05 = new jet::reconstruction_filtering_iir_t(
		0.5F, phenix_acceptance_pseudorapidity_range,
		phenix_acceptance_npixel_pseudorapidity,
		*_background_model_cu_cu_05_perp,
		*_background_model_05_time,
		*_background_model_05_z);
	_kt_03 = new jet::reconstruction_fastjet_t(0.3F);
	_kt_05 = new jet::reconstruction_fastjet_t(0.5F);
	_antikt_03 = new jet::reconstruction_antikt_t(0.3F);
	_antikt_05 = new jet::reconstruction_antikt_t(0.5F);
	_siscone_03 = new jet::reconstruction_siscone_t(0.3F);
	_siscone_05 = new jet::reconstruction_siscone_t(0.5F);
}

event_writer_t::~event_writer_t(void)
{
	if(_output_file != NULL) {
		delete _output_file;
		_output_file = NULL;
	}
	if(_filter_03 != NULL) {
		delete _filter_03;
		_filter_03 = NULL;
	}
	if(_filter_05 != NULL) {
		delete _filter_05;
		_filter_05 = NULL;
	}
	if(_kt_03 != NULL) {
		delete _kt_03;
		_kt_03 = NULL;
	}
	if(_kt_05 != NULL) {
		delete _kt_05;
		_kt_05 = NULL;
	}
	if(_antikt_03 != NULL) {
		delete _antikt_03;
		_antikt_03 = NULL;
	}
	if(_antikt_05 != NULL) {
		delete _antikt_05;
		_antikt_05 = NULL;
	}
	if(_siscone_03 != NULL) {
		delete _siscone_03;
		_siscone_03 = NULL;
	}
	if(_siscone_05 != NULL) {
		delete _siscone_05;
		_siscone_05 = NULL;
	}
	if(_background_model_cu_cu_03_perp != NULL) {
		delete _background_model_cu_cu_03_perp;
		_background_model_cu_cu_03_perp = NULL;
	}
	if(_background_model_cu_cu_05_perp != NULL) {
		delete _background_model_cu_cu_05_perp;
		_background_model_cu_cu_05_perp = NULL;
	}
	if(_background_model_au_au_03_perp != NULL) {
		delete _background_model_au_au_03_perp;
		_background_model_au_au_03_perp = NULL;
	}
	if(_background_model_03_time != NULL) {
		delete _background_model_03_time;
		_background_model_03_time = NULL;
	}
	if(_background_model_03_z != NULL) {
		delete _background_model_03_z;
		_background_model_03_z = NULL;
	}
	if(_background_model_05_time != NULL) {
		delete _background_model_05_time;
		_background_model_05_time = NULL;
	}
	if(_background_model_05_z != NULL) {
		delete _background_model_05_z;
		_background_model_05_z = NULL;
	}
}

int event_writer_t::Init(PHCompositeNode *top_node)
{
	if(!_statistics_output_filename.empty()) {
		_statistics_output_file =
			fopen(_statistics_output_filename.c_str(), "w");
		if(_statistics_output_file == NULL) {
			fprintf(stderr, "%s:%d: WARNING!!!", __FILE__, __LINE__);
			perror("fopen");
		}
	}

	std::vector<jet::histogram_t<float, uint32_t>::fixed_bin_t>
		ert_jet_histogram_bin;

		// Dimension 0
		ert_jet_histogram_bin.push_back(
			jet::histogram_t<float, uint32_t>::fixed_bin_t(
				"centrality (%)",
				nbin_centrality, 1, 101));
		// Dimension 1
		ert_jet_histogram_bin.push_back(
			jet::histogram_t<float, uint32_t>::fixed_bin_t(
				"vertex height (cm)",
				nbin_vertex, -20, 20));
		// Dimension 2
		ert_jet_histogram_bin.push_back(
			jet::histogram_t<float, uint32_t>::fixed_bin_t(
				"pseudorapidity edge distance",
				nbin_edge_distance, 0, 0.15));
		// Dimension 3
		ert_jet_histogram_bin.push_back(
			jet::histogram_t<float, uint32_t>::fixed_bin_t(
				"azimuth edge distance",
				nbin_edge_distance, 0, 0.15));
		// Dimension 4
		ert_jet_histogram_bin.push_back(
			jet::histogram_t<float, uint32_t>::fixed_bin_t(
				"log_10 fake rejection ((GeV/c)^2)",
				nbin_log_fake_rejection, 0.5, 2.0));

		// Dimension 5
		ert_jet_histogram_bin.push_back(
			jet::histogram_t<float, uint32_t>::fixed_bin_t(
				"log_10 jet transverse momentum (GeV/c)",
				nbin_log_perp_or_energy, -1, 2));
		// Dimension 6
		ert_jet_histogram_bin.push_back(
			jet::histogram_t<float, uint32_t>::fixed_bin_t(
				"ERT mode",
				nbin_mode, 0, nbin_mode));
		// Dimension 7
		ert_jet_histogram_bin.push_back(
			jet::histogram_t<float, uint32_t>::fixed_bin_t(
				"supermodule distance",
				nbin_supermodule_distance, 0, 0.5));

		_ert_jet_histogram =
			new jet::histogram_t<float, uint32_t>(
			"jet ERT probability",
			ert_jet_histogram_bin);

		std::vector<jet::histogram_t<float, uint32_t>::fixed_bin_t>
			ert_supermodule_bin;

		// Dimension 0
		ert_supermodule_bin.push_back(
			jet::histogram_t<float, uint32_t>::fixed_bin_t(
				"centrality (%)",
				nbin_centrality, 1, 101));
		// Dimension 1
		ert_supermodule_bin.push_back(
			jet::histogram_t<float, uint32_t>::fixed_bin_t(
				"vertex height (cm)",
				nbin_vertex, -20, 20));

		// Dimension 2
		ert_supermodule_bin.push_back(
			jet::histogram_t<float, uint32_t>::fixed_bin_t(
				"log_10 cluster energy (GeV/c)",
				nbin_log_perp_or_energy, -1, 2));
		// Dimension 3
		ert_supermodule_bin.push_back(
			jet::histogram_t<float, uint32_t>::fixed_bin_t(
				"supermodule",
				nbin_supermodule, 0, nbin_supermodule));
		// Dimension 4
		ert_supermodule_bin.push_back(
			jet::histogram_t<float, uint32_t>::fixed_bin_t(
				"ERT mode",
				nbin_mode, 0, nbin_mode));

		_ert_supermodule_histogram =
			new jet::histogram_t<float, uint32_t>(
			"supermodule ERT probability",
			ert_supermodule_bin);

		_ert_4x4_histogram =
			new jet::histogram_t<float, uint32_t>(
			"supermodule ERT probability",
			ert_supermodule_bin);

		_ert_2x2_histogram =
			new jet::histogram_t<float, uint32_t>(
			"supermodule ERT probability",
			ert_supermodule_bin);

		_ert_supermodule_histogram_cluster =
			new jet::histogram_t<float, uint32_t>(
			"supermodule ERT probability",
			ert_supermodule_bin);

		_ert_4x4_histogram_cluster =
			new jet::histogram_t<float, uint32_t>(
			"supermodule ERT probability",
			ert_supermodule_bin);

		_ert_2x2_histogram_cluster =
			new jet::histogram_t<float, uint32_t>(
			"supermodule ERT probability",
			ert_supermodule_bin);
	
	return 0;
}


int event_writer_t::InitRun(PHCompositeNode *top_node)
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

	void event_writer_t::fill_supermodule(
		const std::vector<jet::ert_hit_t> &ert_hit,
		const jet::event_t &event,
		const std::vector<jet::cluster_t> &cluster)
	{
		float energy_sum_2x2[nbin_2x2];
		float energy_sum_4x4[nbin_2x2];
		float energy_144[nbin_supermodule];
		float energy_2x2[nbin_supermodule];
		float energy_4x4[nbin_supermodule];

		for(unsigned int i = 0; i < nbin_supermodule; i++) {
			energy_144[i] = 0;
			energy_2x2[i] = 0;
			energy_4x4[i] = 0;
		}
		for(unsigned int i = 0; i < nbin_2x2; i++) {
			energy_sum_2x2[i] = 0;
			energy_sum_4x4[i] = 0;
		}

		// 2x2 and 144 tower (supermodule) sums
		for(std::vector<jet::cluster_t>::const_iterator iterator =
				cluster.begin();
			iterator != cluster.end(); iterator++) {
			const jet::tower_id_t tower_id =
				iterator->central_tower_id();
			const int tower_sector = tower_id.sector();
			const int tower_y = tower_id.y();
			const int tower_z = tower_id.z();
			const int index_2x2 = tower_id.lead_scintillator() ?
				(36 * 18) * tower_sector +
				36 * (tower_y >> 1) + (tower_z >> 1) :
				(6 * 36 * 18) + (48 * 24) * (tower_sector - 6) +
				48 * (tower_y >> 1) + (tower_z >> 1);
			energy_sum_2x2[index_2x2] += iterator->energy();
			const int index_supermodule = tower_id.supermodule();
			energy_144[index_supermodule] += iterator->energy();
		}
		// Obtain the 4x4 sum
		for(unsigned int i = 0; i < nbin_2x2; i++) {
			if(i < 6 * 3 * 6 * 36) {
				static const int nz_2x2 = 36;
				static const int ny_2x2 = 18;
				const int z_2x2 = i % ny_2x2;
				const int y_2x2 = (i / ny_2x2) % (nz_2x2 * ny_2x2);
				energy_sum_4x4[i] = energy_sum_2x2[i];
				if(z_2x2 < nz_2x2 - 1) {
					energy_sum_4x4[i] += energy_sum_2x2[i + 1];
				}
				if(y_2x2 < ny_2x2 - 1) {
					energy_sum_4x4[i] += energy_sum_2x2[i + nz_2x2];
				}
				if(z_2x2 < nz_2x2 - 1 && y_2x2 < ny_2x2 - 1) {
					energy_sum_4x4[i] += energy_sum_2x2[i + nz_2x2 + 1];
				}
			}
			else {
				static const int nz_2x2 = 48;
				static const int ny_2x2 = 24;
				const int z_2x2 = i % ny_2x2;
				const int y_2x2 = (i / ny_2x2) % (nz_2x2 * ny_2x2);
				energy_sum_4x4[i] = energy_sum_2x2[i];
				if(z_2x2 < nz_2x2 - 1) {
					energy_sum_4x4[i] += energy_sum_2x2[i + 1];
				}
				if(y_2x2 < ny_2x2 - 1) {
					energy_sum_4x4[i] += energy_sum_2x2[i + nz_2x2];
				}
				if(z_2x2 < nz_2x2 - 1 && y_2x2 < ny_2x2 - 1) {
					energy_sum_4x4[i] += energy_sum_2x2[i + nz_2x2 + 1];
				}
			}
		}
		for(unsigned int i = 0; i < nbin_2x2; i++) {
			const int index_2x2 = i;
			const int sector_2x2 = index_2x2 < 6 * 36 * 18 ?
				index_2x2 / (36 * 18) :
				6 + (index_2x2 - 6 * 36 * 18) / (48 * 24);
			const int y_2x2 = index_2x2 < 6 * 36 * 18 ?
				(index_2x2 / 36) % 18 :
				((index_2x2 - 6 * 36 * 18) / 48) % 24;
			const int z_2x2 = index_2x2 < 6 * 36 * 18 ?
				index_2x2 % 36 :
				(index_2x2 - 6 * 36 * 18) % 48;
			const int index_supermodule = sector_2x2 < 6 ?
				(6 * 3) * sector_2x2 +
				6 * (y_2x2 / 6) + (z_2x2 / 6) :
				(6 * 6 * 3) + (8 * 4) * (sector_2x2 - 6) +
				8 * (y_2x2 / 6) + (z_2x2 / 6);
			energy_4x4[index_supermodule] =
				std::max(energy_4x4[index_supermodule],
						 energy_sum_4x4[i]);
			energy_2x2[index_supermodule] =
				std::max(energy_2x2[index_supermodule],
						 energy_sum_2x2[i]);
		}

		std::vector<float> entry(5, 0);

		entry[0] = event.centrality();
		entry[1] = event.vertex_height().
			beam_beam_counter_mean();

		entry[4] = 0;
		for(unsigned int i = 0; i < nbin_supermodule; i++) {
			entry[3] = i;
			entry[2] = log10f(energy_144[i]);
			reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
				_ert_supermodule_histogram_cluster)->fill(entry, event.scale_down());
			entry[2] = log10f(energy_4x4[i]);
			reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
				_ert_4x4_histogram_cluster)->fill(entry, event.scale_down());
			entry[2] = log10f(energy_2x2[i]);
			reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
				_ert_2x2_histogram_cluster)->fill(entry, event.scale_down());
		}

		for(std::vector<jet::ert_hit_t>::const_iterator iterator =
				ert_hit.begin();
			iterator != ert_hit.end(); iterator++) {
			if(iterator->trigger_mode() !=
			   jet::ert_hit_t::TRIGGER_MODE_RICH) {
				const int id =
					static_cast<int>(iterator->supermodule());

				if(id >= 0 && id < static_cast<int>(nbin_supermodule)) {
					entry[3] = id;
					entry[4] = iterator->trigger_mode() + 1;

					entry[2] = log10f(energy_144[id]);
					reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
						_ert_supermodule_histogram_cluster)->fill(
							entry, event.scale_down());
					entry[2] = log10f(energy_4x4[id]);
					reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
						_ert_4x4_histogram_cluster)->fill(
							entry, event.scale_down());
					entry[2] = log10f(energy_2x2[id]);
					reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
						_ert_2x2_histogram_cluster)->fill(
							entry, event.scale_down());
				}
			}
		}
	}

	void event_writer_t::fill_supermodule_partial(
		const std::vector<jet::ert_hit_t> &ert_hit,
		const jet::event_t &event,
		const std::vector<jet::cluster_t> &cluster)
	{
		float energy_sum_2x2[nbin_2x2];
		float energy_sum_4x4[nbin_2x2];
		float energy_144[nbin_supermodule];
		float energy_2x2[nbin_supermodule];
		float energy_4x4[nbin_supermodule];

		for(unsigned int i = 0; i < nbin_supermodule; i++) {
			energy_144[i] = 0;
			energy_2x2[i] = 0;
			energy_4x4[i] = 0;
		}
		for(unsigned int i = 0; i < nbin_2x2; i++) {
			energy_sum_2x2[i] = 0;
			energy_sum_4x4[i] = 0;
		}

		// 2x2 and 144 tower (supermodule) sums
		for(std::vector<jet::cluster_t>::const_iterator iterator =
				cluster.begin();
			iterator != cluster.end(); iterator++) {
			if(iterator->tower_id().size() !=
			   iterator->partial_energy_sum().size()) {
				std::cerr << __FILE__ << ':' << __LINE__
						  << ": warning: iterator->tower_id().size() "
					"!= iterator->partial_energy_sum().size() ("
						  << iterator->tower_id().size() << " != "
						  << iterator->partial_energy_sum().size()
						  << ')' << std::endl;
			}
			for(unsigned int i = 0; i < iterator->tower_id().size(); i++) {
				const jet::tower_id_t tower_id =
					iterator->tower_id()[i];
				const float energy = i == 0 ?
					iterator->partial_energy_sum()[0] :
					iterator->partial_energy_sum()[i] -
					iterator->partial_energy_sum()[i - 1];
				const int tower_sector = tower_id.sector();
				const int tower_y = tower_id.y();
				const int tower_z = tower_id.z();
				const int index_2x2 = tower_id.lead_scintillator() ?
					(36 * 18) * tower_sector +
					36 * (tower_y >> 1) + (tower_z >> 1) :
					(6 * 36 * 18) + (48 * 24) * (tower_sector - 6) +
					48 * (tower_y >> 1) + (tower_z >> 1);
				energy_sum_2x2[index_2x2] += energy;
				const int index_supermodule = tower_id.supermodule();
				energy_144[index_supermodule] += energy;
			}
		}
		// Obtain the 4x4 sum
		for(unsigned int i = 0; i < nbin_2x2; i++) {
			if(i < 6 * 3 * 6 * 36) {
				static const int nz_2x2 = 36;
				static const int ny_2x2 = 18;
				const int z_2x2 = i % ny_2x2;
				const int y_2x2 = (i / ny_2x2) % (nz_2x2 * ny_2x2);
				energy_sum_4x4[i] = energy_sum_2x2[i];
				if(z_2x2 < nz_2x2 - 1) {
					energy_sum_4x4[i] += energy_sum_2x2[i + 1];
				}
				if(y_2x2 < ny_2x2 - 1) {
					energy_sum_4x4[i] += energy_sum_2x2[i + nz_2x2];
				}
				if(z_2x2 < nz_2x2 - 1 && y_2x2 < ny_2x2 - 1) {
					energy_sum_4x4[i] += energy_sum_2x2[i + nz_2x2 + 1];
				}
			}
			else {
				static const int nz_2x2 = 48;
				static const int ny_2x2 = 24;
				const int z_2x2 = i % ny_2x2;
				const int y_2x2 = (i / ny_2x2) % (nz_2x2 * ny_2x2);
				energy_sum_4x4[i] = energy_sum_2x2[i];
				if(z_2x2 < nz_2x2 - 1) {
					energy_sum_4x4[i] += energy_sum_2x2[i + 1];
				}
				if(y_2x2 < ny_2x2 - 1) {
					energy_sum_4x4[i] += energy_sum_2x2[i + nz_2x2];
				}
				if(z_2x2 < nz_2x2 - 1 && y_2x2 < ny_2x2 - 1) {
					energy_sum_4x4[i] += energy_sum_2x2[i + nz_2x2 + 1];
				}
			}
		}
		for(unsigned int i = 0; i < nbin_2x2; i++) {
			const int index_2x2 = i;
			const int sector_2x2 = index_2x2 < 6 * 36 * 18 ?
				index_2x2 / (36 * 18) :
				6 + (index_2x2 - 6 * 36 * 18) / (48 * 24);
			const int y_2x2 = index_2x2 < 6 * 36 * 18 ?
				(index_2x2 / 36) % 18 :
				((index_2x2 - 6 * 36 * 18) / 48) % 24;
			const int z_2x2 = index_2x2 < 6 * 36 * 18 ?
				index_2x2 % 36 :
				(index_2x2 - 6 * 36 * 18) % 48;
			const int index_supermodule = sector_2x2 < 6 ?
				(6 * 3) * sector_2x2 +
				6 * (y_2x2 / 6) + (z_2x2 / 6) :
				(6 * 6 * 3) + (8 * 4) * (sector_2x2 - 6) +
				8 * (y_2x2 / 6) + (z_2x2 / 6);
			energy_4x4[index_supermodule] =
				std::max(energy_4x4[index_supermodule],
						 energy_sum_4x4[i]);
			energy_2x2[index_supermodule] =
				std::max(energy_2x2[index_supermodule],
						 energy_sum_2x2[i]);
		}

		std::vector<float> entry(5, 0);

		entry[0] = event.centrality();
		entry[1] = event.vertex_height().
			beam_beam_counter_mean();

		entry[4] = 0;
		for(unsigned int i = 0; i < nbin_supermodule; i++) {
			entry[3] = i;
			entry[2] = log10f(energy_144[i]);
			reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
				_ert_supermodule_histogram)->fill(entry, event.scale_down());
			entry[2] = log10f(energy_4x4[i]);
			reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
				_ert_4x4_histogram)->fill(entry, event.scale_down());
			entry[2] = log10f(energy_2x2[i]);
			reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
				_ert_2x2_histogram)->fill(entry, event.scale_down());
		}

		for(std::vector<jet::ert_hit_t>::const_iterator iterator =
				ert_hit.begin();
			iterator != ert_hit.end(); iterator++) {
			if(iterator->trigger_mode() !=
			   jet::ert_hit_t::TRIGGER_MODE_RICH) {
				const int id =
					static_cast<int>(iterator->supermodule());

				if(id >= 0 && id < static_cast<int>(nbin_supermodule)) {
					entry[3] = id;
					entry[4] = iterator->trigger_mode() + 1;

					entry[2] = log10f(energy_144[id]);
					reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
						_ert_supermodule_histogram)->fill(
							entry, event.scale_down());
					entry[2] = log10f(energy_4x4[id]);
					reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
						_ert_4x4_histogram)->fill(
							entry, event.scale_down());
					entry[2] = log10f(energy_2x2[id]);
					reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
						_ert_2x2_histogram)->fill(
							entry, event.scale_down());
				}
			}
		}
	}

	void event_writer_t::fill_jet(
		const jet::jet_t &jet,
		const std::vector<jet::ert_hit_t> &ert_hit,
		const jet::event_t &event)
	{
		std::vector<bool> mode(nbin_mode - 1, false);
		std::vector<float> distance(nbin_mode - 1, NAN);

		for(unsigned int i = 0; i < nbin_mode - 1; i++) {
			distance[i] = NAN;
		}

		for(std::vector<jet::ert_hit_t>::const_iterator iterator =
				ert_hit.begin();
			iterator != ert_hit.end(); iterator++) {
			mode[iterator->trigger_mode()] = true;
			mode[nbin_mode - 2] = true;

			const float distance_pseudorapidity =
				jet.momentum().pseudorapidity() -
				iterator->supermodule().approximate_pseudorapidity();
			const float distance_azimuth =
				jet::angular_range_reduce(
					jet.momentum().azimuth() -
					iterator->supermodule().approximate_azimuth());
			const float distance_radial = sqrtf(
				distance_pseudorapidity * distance_pseudorapidity +
				distance_azimuth * distance_azimuth);

			distance[iterator->trigger_mode()] =
				distance_radial > distance[iterator->trigger_mode()] ?
				distance[iterator->trigger_mode()] : distance_radial;
			distance[nbin_mode - 2] =
				distance_radial > distance[nbin_mode - 2] ?
				distance[nbin_mode - 2] : distance_radial;
		}
		// fprintf(stderr, "%u %f\n", iterator->trigger_mode(), distance[iterator->trigger_mode()]);
		// fprintf(stderr, "%f: %f %f %f %f\n", jet.momentum().perp(), distance[0], distance[1], distance[2], distance[3]);

		std::vector<float> entry(8, 0);

		entry[0] = event.centrality();
		entry[1] = event.vertex_height().
			beam_beam_counter_mean();
		entry[2] =
			jet.momentum().
			pseudorapidity_edge_distance(
				pseudorapidity_min,
				pseudorapidity_max);
		entry[3] =
			jet.momentum().
			azimuth_edge_distance(
				azimuth_west_min, azimuth_west_max,
				azimuth_east_min, azimuth_east_max);
		entry[4] = FLT_MAX;
		entry[5] = log10f(
			jet.momentum().perp());

		entry[6] = 0;
		entry[7] = 0;
		reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
			_ert_jet_histogram)->fill(entry, event.scale_down());
		for(unsigned int i = 0; i < nbin_mode; i++) {
			if(mode[i]) {
				entry[6] = i + 1;
				entry[7] = distance[i];
				reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
					_ert_jet_histogram)->fill(
						entry, event.scale_down());
			}
		}
	}

int event_writer_t::process_event(PHCompositeNode *top_node)
{
	const EventHeader *phenix_header =
		findNode::getClass<EventHeader>(top_node, "EventHeader");

	// No header is considered fatal
	if(phenix_header == NULL) {
		std::cerr << __FILE__ << ':' << __LINE__
				  << ": phenix_header = NULL" << std::endl;
		std::cerr << ">>> End: " << __FUNCTION__ << std::endl;
		return 1;
	}

	const TrigLvl1 *phenix_level_1_trigger =
		findNode::getClass<TrigLvl1>(top_node, "TrigLvl1");
	uint64_t time = 0U;
	uint32_t gl1_raw = 0U;
	uint32_t gl1_live = 0U;
	uint32_t gl1_scaled = 0U;

	if(phenix_level_1_trigger != NULL) {
		time = static_cast<uint64_t>(
				phenix_level_1_trigger->get_lvl1_beam_clk(0)) |
			(static_cast<uint64_t>(
				phenix_level_1_trigger->get_lvl1_beam_clk(1)) << 32);
		gl1_raw = phenix_level_1_trigger->get_lvl1_trigraw();
		gl1_live = phenix_level_1_trigger->get_lvl1_triglive();
		gl1_scaled = phenix_level_1_trigger->get_lvl1_trigscaled();
	}

	const jet::gl1_state_t gl1_state =
		jet::gl1_state_t(gl1_raw, gl1_live, gl1_scaled);

	const ErtOut *phenix_ert =
		findNode::getClass<ErtOut>(top_node, "ErtOut");
	const PHGlobal *phenix_global =
		findNode::getClass<PHGlobal>(top_node, "PHGlobal");

	jet::vertex_height_t vertex_height(NAN, NAN, NAN, NAN);
	float centrality = NAN;
	float bbc_t0 = NAN;

	if(phenix_global != NULL) {
		vertex_height = jet::vertex_height_t(
			phenix_global->getBbcZVertex(),
			phenix_global->getBbcZVertexError(),
			phenix_global->getZdcZVertex(),
			phenix_global->getZdcZVertexError());
		centrality = phenix_global->getCentrality();
		bbc_t0 = phenix_global->getBbcTimeZero();
	}

	const uint64_t event_number = phenix_header->get_EvtSequence();
	std::vector<jet::ert_hit_t> ert_hit;

	if(phenix_ert != NULL) {
		const unsigned int nert_hit = phenix_ert->get_ERThit_N();

		ert_hit.reserve(nert_hit);
		for(unsigned int i = 0; i < nert_hit; i++) {
			const unsigned int mode1 = phenix_ert->get_ERTtrigmode(i);
			const unsigned int arm1 = phenix_ert->get_ERTarm(i);
			const unsigned int sector1 = phenix_ert->get_ERTsector(i);
			const unsigned int supermodule1 = phenix_ert->get_ERTsm(i);

			ert_hit.push_back(jet::ert_hit_t(
				mode1,
				jet::supermodule_id_t(
					arm1, sector1, supermodule1, mode1 == 4)));
		}
	}

	const ReactionPlaneObject *phenix_reaction_plane =
		findNode::getClass<ReactionPlaneObject>(
			top_node, "ReactionPlaneObject");
	const RpSumXYObject *phenix_reaction_plane_flow_q =
		findNode::getClass<RpSumXYObject>(
			top_node, "RpSumXYObject");

	jet::phenix_reaction_plane_t reaction_plane;

	// Reaction plane is optional
	if(phenix_reaction_plane != NULL &&
	   phenix_reaction_plane_flow_q != NULL) {
		if(_legacy) {
			reaction_plane = jet::phenix_reaction_plane_t(
				// BBC
				jet::north_south_reaction_plane_t(
					jet::reaction_plane_t(
						phenix_reaction_plane->getBBCrp12(),
						jet::cartesian_point_t(
							phenix_reaction_plane_flow_q->getBBCsumX12(),
							phenix_reaction_plane_flow_q->getBBCsumY12()),
						phenix_reaction_plane_flow_q->getBBCsumW2()),
					jet::reaction_plane_t(
						phenix_reaction_plane->getBBCrp11(),
						jet::cartesian_point_t(
							phenix_reaction_plane_flow_q->getBBCsumX11(),
							phenix_reaction_plane_flow_q->getBBCsumY11()),
						phenix_reaction_plane_flow_q->getBBCsumW1()),
					jet::reaction_plane_t(
						phenix_reaction_plane->getBBCrp10(),
						jet::cartesian_point_t(
							phenix_reaction_plane_flow_q->getBBCsumX10(),
							phenix_reaction_plane_flow_q->getBBCsumY10()),
						phenix_reaction_plane_flow_q->getBBCsumW0())),
				jet::north_south_reaction_plane_t(),
				jet::north_south_reaction_plane_t(),
				jet::north_south_reaction_plane_t(),
				jet::north_south_reaction_plane_t());
		}
		else {
			reaction_plane = jet::phenix_reaction_plane_t(
				// BBC
				jet::north_south_reaction_plane_t(
					jet::reaction_plane_t(
						phenix_reaction_plane->getBBCrp12(),
						jet::cartesian_point_t(
							phenix_reaction_plane_flow_q->getBBCsumX12(),
							phenix_reaction_plane_flow_q->getBBCsumY12()),
						phenix_reaction_plane_flow_q->getBBCsumW2()),
					jet::reaction_plane_t(
						phenix_reaction_plane->getBBCrp11(),
						jet::cartesian_point_t(
							phenix_reaction_plane_flow_q->getBBCsumX11(),
							phenix_reaction_plane_flow_q->getBBCsumY11()),
						phenix_reaction_plane_flow_q->getBBCsumW1()),
					jet::reaction_plane_t(
						phenix_reaction_plane->getBBCrp10(),
						jet::cartesian_point_t(
							phenix_reaction_plane_flow_q->getBBCsumX10(),
							phenix_reaction_plane_flow_q->getBBCsumY10()),
						phenix_reaction_plane_flow_q->getBBCsumW0())),
				// RXNP combined
				jet::north_south_reaction_plane_t(
					jet::reaction_plane_t(
						phenix_reaction_plane->getRXNrp18(),
						jet::cartesian_point_t(
							phenix_reaction_plane_flow_q->getRXNsumX18(),
							phenix_reaction_plane_flow_q->getRXNsumY18()),
						phenix_reaction_plane_flow_q->getRXNsumW8()),
					jet::reaction_plane_t(
						phenix_reaction_plane->getRXNrp15(),
						jet::cartesian_point_t(
							phenix_reaction_plane_flow_q->getRXNsumX15(),
							phenix_reaction_plane_flow_q->getRXNsumY15()),
						phenix_reaction_plane_flow_q->getRXNsumW5()),
					jet::reaction_plane_t(
						phenix_reaction_plane->getRXNrp12(),
						jet::cartesian_point_t(
							phenix_reaction_plane_flow_q->getRXNsumX12(),
							phenix_reaction_plane_flow_q->getRXNsumY12()),
						phenix_reaction_plane_flow_q->getRXNsumW2())),
				// RXNP inner
				jet::north_south_reaction_plane_t(
					jet::reaction_plane_t(
						phenix_reaction_plane->getRXNrp16(),
						jet::cartesian_point_t(
							phenix_reaction_plane_flow_q->getRXNsumX16(),
							phenix_reaction_plane_flow_q->getRXNsumY16()),
						phenix_reaction_plane_flow_q->getRXNsumW6()),
					jet::reaction_plane_t(
						phenix_reaction_plane->getRXNrp13(),
						jet::cartesian_point_t(
							phenix_reaction_plane_flow_q->getRXNsumX13(),
							phenix_reaction_plane_flow_q->getRXNsumY13()),
						phenix_reaction_plane_flow_q->getRXNsumW3()),
					jet::reaction_plane_t(
						phenix_reaction_plane->getRXNrp10(),
						jet::cartesian_point_t(
							phenix_reaction_plane_flow_q->getRXNsumX10(),
							phenix_reaction_plane_flow_q->getRXNsumY10()),
						phenix_reaction_plane_flow_q->getRXNsumW0())),
				// RXNP outer
				jet::north_south_reaction_plane_t(
					jet::reaction_plane_t(
						phenix_reaction_plane->getRXNrp17(),
						jet::cartesian_point_t(
							phenix_reaction_plane_flow_q->getRXNsumX17(),
							phenix_reaction_plane_flow_q->getRXNsumY17()),
						phenix_reaction_plane_flow_q->getRXNsumW7()),
					jet::reaction_plane_t(
						phenix_reaction_plane->getRXNrp14(),
						jet::cartesian_point_t(
							phenix_reaction_plane_flow_q->getRXNsumX14(),
							phenix_reaction_plane_flow_q->getRXNsumY14()),
						phenix_reaction_plane_flow_q->getRXNsumW4()),
					jet::reaction_plane_t(
						phenix_reaction_plane->getRXNrp11(),
						jet::cartesian_point_t(
							phenix_reaction_plane_flow_q->getRXNsumX11(),
							phenix_reaction_plane_flow_q->getRXNsumY11()),
						phenix_reaction_plane_flow_q->getRXNsumW1())),
				// MPC
				jet::north_south_reaction_plane_t(
					jet::reaction_plane_t(
						phenix_reaction_plane->getMPCrp12(),
						jet::cartesian_point_t(
							phenix_reaction_plane_flow_q->getMPCsumX12(),
							phenix_reaction_plane_flow_q->getMPCsumY12()),
						phenix_reaction_plane_flow_q->getMPCsumW2()),
					jet::reaction_plane_t(
						phenix_reaction_plane->getMPCrp11(),
						jet::cartesian_point_t(
							phenix_reaction_plane_flow_q->getMPCsumX11(),
							phenix_reaction_plane_flow_q->getMPCsumY11()),
						phenix_reaction_plane_flow_q->getMPCsumW1()),
					jet::reaction_plane_t(
						phenix_reaction_plane->getMPCrp10(),
						jet::cartesian_point_t(
							phenix_reaction_plane_flow_q->getMPCsumX10(),
							phenix_reaction_plane_flow_q->getMPCsumY10()),
						phenix_reaction_plane_flow_q->getMPCsumW0())));
		}
	}

	const PHCentralTrack *phenix_central =
		findNode::getClass<PHCentralTrack>(
			top_node, "PHCentralTrack");
	std::vector<jet::track_t> track;

	if(phenix_central == NULL && _nevent < 20) {
		std::cerr << __FILE__ << ':' << __LINE__ << ": "
			"information: phenix_central == NULL" << std::endl;
	}
	else if(phenix_central != NULL &&
			phenix_central->get_npart() == 0 && _nevent < 20) {
		std::cerr << __FILE__ << ':' << __LINE__ << ": "
			"information: phenix_central->get_npart() == 0" << std::endl;
	}

	if(phenix_central != NULL) {
		for(int track_idx = 0;
			track_idx < (int)phenix_central->get_npart(); track_idx++) {
			const int quality =
				phenix_central->get_quality(track_idx);

			// Drop tracks by quality
			if(quality <= 7) {
				continue;
			}

			// Momentum arithmetics
			const double polar_angle =
				phenix_central->get_the0(track_idx);

			// Drop illegal track
			if(polar_angle <= -9999) {
				continue;
			}

			I(polar_angle <= -9999 ||
			  (polar_angle >= 0 && polar_angle <= M_PI));

			const double momentum_norm =
				phenix_central->get_mom(track_idx);

			I(momentum_norm >= 0);

			const float perp =
				(float)(momentum_norm * sin(polar_angle));

			I(perp >= 0);

			const float pseudorapidity =
				(float)(-log(tan(0.5 * polar_angle)));
			const float azimuth =
				phenix_central->get_phi0(track_idx);

			// Drop illegal track
			if(azimuth <= -9999) {
				continue;
			}

			I(azimuth >= -M_PI && azimuth <= 2 * M_PI);

			const float azimuth_reduced =
				jet::angular_range_reduce(azimuth);

			I(azimuth_reduced > -M_PI &&
			  azimuth_reduced <= M_PI);

			// Other information

			// PC3 first, since tracks not satisfying PC3 matching cuts
			// are thrown out immediately
			const float pc3_sigma_displacement_azimuth =
				phenix_central->get_pc3sdphi(track_idx);
			const float pc3_sigma_displacement_height =
				phenix_central->get_pc3sdz(track_idx);

			// Apply PC3 Matching cut

			if(pc3_sigma_displacement_azimuth <= -9999 ||
			   pc3_sigma_displacement_height <= -9999 ||
			   !(fabs(pc3_sigma_displacement_azimuth) <= 3 &&
				 fabs(pc3_sigma_displacement_height) <= 3)) {
				continue;
			}

			const int charge =
				phenix_central->get_charge(track_idx);
			const int cluster_id =
				phenix_central->get_emcid(track_idx);
			const float dc_crossing_polar_angle =
				phenix_central->get_beta(track_idx);
			const float dc_crossing_azimuth =
				phenix_central->get_phi(track_idx);

			I(dc_crossing_polar_angle >= 0 &&
			  dc_crossing_polar_angle <= M_PI);
			I(dc_crossing_azimuth > -M_PI &&
			  dc_crossing_azimuth <= 2 * M_PI);

			const float pc1_crossing_height =
				phenix_central->get_zed(track_idx);
			const int rich_ring_nnormal_area_phototube =
				phenix_central->get_n0(track_idx);
			const float rich_ring_nnormal_area_photoelectron =
				phenix_central->get_npe0(track_idx);
			const int rich_ring_nwide_area_phototube =
				phenix_central->get_n1(track_idx);
			const float rich_ring_nwide_area_photoelectron =
				phenix_central->get_npe1(track_idx);
			const float rich_ring_unitful_chi_square =
				phenix_central->get_chi2(track_idx);
			const float rich_ring_displacement =
				phenix_central->get_disp(track_idx);

			I(F(rich_ring_unitful_chi_square));
			I(F(rich_ring_displacement));

			const float pc2_sigma_displacement_azimuth =
				phenix_central->get_pc2sdphi(track_idx);
			const float pc2_sigma_displacement_height =
				phenix_central->get_pc2sdz(track_idx);
			const float tec_sigma_displacement_azimuth =
				phenix_central->get_tecsdphi(track_idx);
			const float tec_sigma_displacement_inclination =
				phenix_central->get_tecsdalpha(track_idx);
			const float cluster_sigma_displacement_azimuth =
				phenix_central->get_emcsdphi(track_idx);
			const float cluster_sigma_displacement_height =
				phenix_central->get_emcsdz(track_idx);

			const float dc_crossing_azimuth_reduced =
				dc_crossing_azimuth == -9999 ? dc_crossing_azimuth :
				jet::angular_range_reduce(dc_crossing_azimuth);

			I(dc_crossing_azimuth_reduced > -M_PI &&
			  dc_crossing_azimuth_reduced <= M_PI);

			jet::track_t track1(
				NAN, perp, pseudorapidity, azimuth_reduced, charge,
				quality, cluster_id,
				dc_crossing_polar_angle,
				dc_crossing_azimuth_reduced,
				pc1_crossing_height,
				rich_ring_nnormal_area_phototube,
				rich_ring_nnormal_area_photoelectron,
				rich_ring_nwide_area_phototube,
				rich_ring_nwide_area_photoelectron,
				rich_ring_unitful_chi_square,
				rich_ring_displacement,
				pc2_sigma_displacement_azimuth,
				pc2_sigma_displacement_height,
				tec_sigma_displacement_azimuth,
				tec_sigma_displacement_inclination,
				pc3_sigma_displacement_azimuth,
				pc3_sigma_displacement_height,
				cluster_sigma_displacement_azimuth,
				cluster_sigma_displacement_height);

			track.push_back(track1);
		}
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

	bool write_event;
	unsigned int scale_down = UINT_MAX;
	static const unsigned int nfake_rejection_level = 3;
	std::vector<jet::jet_t> jet_03;
	std::vector<jet::jet_t> jet_05;
	std::vector<jet::jet_t> jet_kt_03;
	std::vector<jet::jet_t> jet_kt_05;
	std::vector<jet::jet_t> jet_antikt_03;
	std::vector<jet::jet_t> jet_antikt_05;
	std::vector<jet::jet_t> jet_siscone_03;
	std::vector<jet::jet_t> jet_siscone_05;

	if(_mode == SCALEDOWN_NONE) {
		write_event = true;
		scale_down = 1;
	}
	else if(_mode == SCALEDOWN_FIXED_4) {
		write_event = _nevent % 4 == 0;
		scale_down = 4;
	}
	else {
		jet::event_t test_event(
			event_number, time, 1.0F, track, cluster,
			std::vector<jet::cluster_t>(),
			std::vector<jet::particle_t>(),
			ert_hit, gl1_state, vertex_height, bbc_t0,
			centrality, reaction_plane);

		// Test event
		bool bad_event;
		float run_timing_mean[8];
		float tower_timing_mean[jet::tower_id_t::ntower];
		float tower_timing_standard_deviation[jet::tower_id_t::ntower];

		switch(_mode) {
		case JET_RECONSTRUCT_RUN_5_P_P:
#include <extract/run_5_p_p_tower_map.h>
#include <extract/run_5_p_p_tower_scale.h>
			test_event = test_event.phenix_reconstruct(
				bad_event,
				run_5_p_p_tower_map_3x3, run_5_p_p_tower_scale,
				NULL, NULL, NULL,
				false, true, false, false);
			break;
		case JET_RECONSTRUCT_RUN_5_CU_CU:
#include <extract/run_5_cu_cu_tower_map.h>
#include <extract/run_5_cu_cu_tower_scale.h>
			test_event = test_event.phenix_reconstruct(
				bad_event,
				run_5_cu_cu_tower_map_3x3, run_5_cu_cu_tower_scale,
				NULL, NULL, NULL,
				false, true, false, false);
			break;
		case JET_RECONSTRUCT_RUN_7_AU_AU:
#include <extract/run_7_au_au_tower_map.h>
#include <extract/run_7_au_au_tower_scale.h>
			test_event = test_event.phenix_reconstruct(
				bad_event,
				run_7_au_au_tower_map_3x3, run_7_au_au_tower_scale,
				NULL, NULL, NULL,
				false, true, false, false);
			break;
		case JET_RECONSTRUCT_RUN_8_D_AU:
#include <extract/run_8_d_au_tower_map.h>
#include <extract/run_8_d_au_tower_scale.h>
			test_event = test_event.phenix_reconstruct(
				bad_event,
				run_5_cu_cu_tower_map_3x3, run_5_cu_cu_tower_scale,
				NULL, NULL, NULL,
				false, true, false, false);
			break;
		}

		if(_mode == JET_RECONSTRUCT_RUN_5_P_P) {
			test_event.centrality() = NAN;
		}

		const jet::collision_geometry_t geometry(
			vertex_height.beam_beam_counter_mean(),
			test_event.centrality(),
			test_event.reaction_plane());

		jet_03 = _filter_03->reconstruct(
			test_event.track(), geometry);
		if(_mode != JET_RECONSTRUCT_RUN_7_AU_AU) {
			jet_05 = _filter_05->reconstruct(
				test_event.track(), geometry);
		}
		jet_kt_03 = _kt_03->reconstruct(
			test_event.track(), geometry);
		jet_kt_05 = _kt_05->reconstruct(
			test_event.track(), geometry);
		jet_antikt_03 = _antikt_03->reconstruct(
			test_event.track(), geometry);
		jet_antikt_05 = _antikt_05->reconstruct(
			test_event.track(), geometry);
		jet_siscone_03 = _siscone_03->reconstruct(
			test_event.track(), geometry);
		jet_siscone_05 = _siscone_05->reconstruct(
			test_event.track(), geometry);

		static const float quadratic_fake_rejection_level_cu_cu[
			nfake_rejection_level] = {
			4.0, 8.0, 16.0
		};
		static const float quadratic_fake_rejection_level_au_au[
			nfake_rejection_level] = {
			8.0, 16.0, 32.0
		};
		const float *plevel = NULL;

		switch(_mode) {
		case JET_RECONSTRUCT_RUN_5_CU_CU:
			plevel = quadratic_fake_rejection_level_cu_cu;
			break;
		case JET_RECONSTRUCT_RUN_7_AU_AU:
			plevel = quadratic_fake_rejection_level_au_au;
			break;
		}

		std::vector<jet::jet_t> jet_03_fr[nfake_rejection_level];
		std::vector<jet::jet_t> jet_05_fr[nfake_rejection_level];

		if(plevel != NULL) {
			for(unsigned int i = 0; i < nfake_rejection_level; i++) {
				jet_03_fr[i] = gaussian_fake_rejection(
					jet_03, test_event, plevel[i]);
				jet_05_fr[i] = gaussian_fake_rejection(
					jet_05, test_event, plevel[i]);
			}
		}

		static const unsigned int scale_down_max_p_p = 32U;
		static const unsigned int scale_down_max_cu_cu = 64U;
		static const unsigned int scale_down_max_au_au = 64U;
		static const unsigned int nthreshold_p_p = 3;
		static const float perp_threshold_p_p[nthreshold_p_p] = {
			4.0F, 1e+6F, 1e+6F
		};
		static const unsigned int nthreshold_cu_cu = 4;
		static const float perp_threshold_cu_cu[
			nthreshold_cu_cu * nfake_rejection_level] = {
			8.0F, 1e+6F, 1e+6F, 1e+6F,
			4.0F, 1e+6F, 1e+6F, 1e+6F,
			4.0F, 1e+6F, 1e+6F, 1e+6F
		};
		static const unsigned int nthreshold_au_au = 4;
		static const float perp_threshold_au_au[
			nthreshold_au_au * nfake_rejection_level] = {
			64.0F, 32.0F, 16.0F, 8.0F,
			32.0F, 16.0F, 8.0F, 4.0F,
			16.0F, 8.0F, 4.0F, 2.0F
		};
		unsigned int nthreshold;
		const float *pthreshold;

		switch(_mode) {
		case JET_RECONSTRUCT_RUN_5_P_P:
		case JET_RECONSTRUCT_RUN_8_D_AU:
			scale_down = scale_down_max_p_p;
			nthreshold = nthreshold_p_p;
			pthreshold = perp_threshold_p_p;
			if(!jet_03.empty()) {
				float max_perp = 0;

				for(std::vector<jet::jet_t>::const_iterator
						iterator = jet_03.begin();
					iterator != jet_03.end(); iterator++) {
					max_perp = std::max(
						max_perp, iterator->momentum().perp());
				}
				for(unsigned int i = 0; i < nthreshold; i++) {
					if(max_perp >= pthreshold[i]) {
						scale_down = std::min(scale_down, 1U << i);
						break;
					}
				}
			}
			break;
		case JET_RECONSTRUCT_RUN_5_CU_CU:
			scale_down = scale_down_max_cu_cu;
			nthreshold = nthreshold_cu_cu;
			pthreshold = perp_threshold_cu_cu;
			if(!jet_03.empty()) {
				float max_perp = 0;

				for(unsigned int j = 0; j < nfake_rejection_level; j++) {
					for(std::vector<jet::jet_t>::const_iterator
							iterator = jet_03_fr[j].begin();
						iterator != jet_03_fr[j].end(); iterator++) {
						max_perp = std::max(
							max_perp, iterator->momentum().perp());
					}
					for(unsigned int i = 0; i < nthreshold; i++) {
						if(max_perp >= pthreshold[
								j * nthreshold + i]) {
							scale_down = std::min(scale_down, 1U << i);
							break;
						}
					}
				}
			}
			break;
		case JET_RECONSTRUCT_RUN_7_AU_AU:
			scale_down = scale_down_max_au_au;
			nthreshold = nthreshold_au_au;
			pthreshold = perp_threshold_au_au;
			if(!jet_03.empty()) {
				float max_perp = 0;

				for(unsigned int j = 0; j < nfake_rejection_level; j++) {
					for(std::vector<jet::jet_t>::const_iterator
							iterator = jet_03_fr[j].begin();
						iterator != jet_03_fr[j].end(); iterator++) {
						max_perp = std::max(
							max_perp, iterator->momentum().perp());
					}
					for(unsigned int i = 0; i < nthreshold; i++) {
						if(max_perp >= pthreshold[
								j * nthreshold + i]) {
							scale_down = std::min(scale_down, 1U << i);
							break;
						}
					}
				}
			}
			break;
		}

		write_event = (_nevent % scale_down) == 0;
	}

	jet::event_t event(
		event_number, time, scale_down, track, cluster,
		g_cluster_raw,
		std::vector<jet::particle_t>(),
		ert_hit, gl1_state, vertex_height, bbc_t0,
		centrality, reaction_plane);

	if(write_event) {
		if(_nevent % (_nevent >= 1000 ? 1000 : _nevent >= 100 ? 100 :
					  _nevent >= 20 ? 10 : 1) == 0) {
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": information: experiment_run_number = "
					  << _experiment_run_number << ", nevent = "
					  << _nevent << ", centrality = " << centrality
					  << ", jet multiplicity: sigma = 0.3: "
					  << jet_03.size() << ", sigma = 0.5: "
					  << jet_05.size() << ", d = 0.3 (kT): "
					  << jet_kt_03.size() << ", d = 0.5 (kT): "
					  << jet_kt_05.size() << ", d = 0.3 (anti-kT): "
					  << jet_antikt_03.size() << ", d = 0.5 (anti-kT): "
					  << jet_antikt_05.size() << ", R = 0.3: "
					  << jet_siscone_03.size() << ", R = 0.5: "
					  << jet_siscone_05.size() << ", scale_down = "
					  << scale_down << std::endl;
		}

		event.parton().clear();
		event.parton().reserve(jet_03.size());
		for(std::vector<jet::jet_t>::const_iterator iterator =
				jet_03.begin();
			iterator != jet_03.end(); iterator++) {
			event.parton().push_back(
				jet::particle_t(iterator->momentum()));
		}

		*_output_file << event;
	}

	event.scale_down() = 1;
	if(!event.cluster().empty()) {
		fill_supermodule(event.ert_hit(), event, event.cluster());
		fill_supermodule_partial(event.ert_hit(), event, event.cluster());
	}
	if(!jet_03.empty()) {
		fill_jet(jet_03[0], event.ert_hit(), event);
	}

	_nevent++;

	if(_nevent < 20) {
		std::cerr << __FILE__ << ':' << __LINE__ << ": " << std::endl;
	}

	return 0;
}

int event_writer_t::End(PHCompositeNode *top_node)
{
	fprintf(stderr, "%s:%d: %s %lu\n", __FILE__, __LINE__,
			_statistics_output_filename.c_str(),
			(unsigned long)_statistics_output_file);
	if(_statistics_output_file != NULL) {
		jet::fwrite_histogram(
			*reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
				_ert_jet_histogram),
			_statistics_output_file);
		jet::fwrite_histogram(
			*reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
				_ert_supermodule_histogram),
			_statistics_output_file);
		jet::fwrite_histogram(
			*reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
				_ert_4x4_histogram),
			_statistics_output_file);
		jet::fwrite_histogram(
			*reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
				_ert_2x2_histogram),
			_statistics_output_file);
		jet::fwrite_histogram(
			*reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
				_ert_supermodule_histogram_cluster),
			_statistics_output_file);
		jet::fwrite_histogram(
			*reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
				_ert_4x4_histogram_cluster),
			_statistics_output_file);
		jet::fwrite_histogram(
			*reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
				_ert_2x2_histogram_cluster),
			_statistics_output_file);
		fclose(_statistics_output_file);
	}
	delete reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
		_ert_jet_histogram);
	delete reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
		_ert_supermodule_histogram);
	delete reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
		_ert_4x4_histogram);
	delete reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
		_ert_2x2_histogram);
	delete reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
		_ert_supermodule_histogram_cluster);
	delete reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
		_ert_4x4_histogram_cluster);
	delete reinterpret_cast<jet::histogram_t<float, uint32_t> *>(
		_ert_2x2_histogram_cluster);
	return 0;
}
