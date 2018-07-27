#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H
#include <stdint.h>
#include <jetevent/event.h>

namespace jet {

	const unsigned int supermodule_id_t::
	_lead_scintillator_sector_stride =
		phenix_calorimetry_t::_lead_scintillator_sector_stride / (12 * 12);

	const unsigned int supermodule_id_t::
	_lead_glass_sector_stride =
		phenix_calorimetry_t::_lead_glass_sector_stride / (12 * 12);

	const unsigned int supermodule_id_t::
	_lead_scintillator_z_stride =
		phenix_calorimetry_t::_lead_scintillator_z_stride / 12;

	const unsigned int supermodule_id_t::
	_lead_glass_z_stride =
		phenix_calorimetry_t::_lead_glass_z_stride / 12;

	const unsigned int supermodule_id_t::lead_scintillator_offset =
		0;

	const unsigned int supermodule_id_t::lead_glass_offset =
		_nsector_lead_scintillator *
		_lead_scintillator_sector_stride;

	const unsigned int supermodule_id_t::rich_offset =
		_nsector_lead_scintillator *
		_lead_scintillator_sector_stride +
		_nsector_lead_glass * _lead_glass_sector_stride;

	const unsigned int supermodule_id_t::nsupermodule =
		_nsector_lead_scintillator *
		_lead_scintillator_sector_stride +
		_nsector_lead_glass * _lead_glass_sector_stride +
		_nsector_rich * _rich_sector_stride;

	std::string gl1_state_t::binary_str(const uint32_t n) const
	{
		std::string ret(32, '0');

		for(int i = 0; i < 32; i++) {
			if((n & (1U << i)) != 0) {
				ret[31 - i] = '1';
			}
		}

		return ret;
	}

	void supermodule_id_t::
	set_supermodule_id(const int arm, const int sector,
					   const int reduced_id, const bool rich)
	{
		I(arm >= 0 && arm < 2);
		I(sector >= 0 && sector < 4);

		const int continuous_sector = arm == 0 ? sector :
			sector < 2 ? sector + 6 : sector + 2;

		I(continuous_sector >= 0 && continuous_sector < 8);

		I(lead_scintillator_offset == 0);
		I(_lead_scintillator_sector_stride == 18);
		I(lead_glass_offset == 6 * 18);
		I(_lead_glass_sector_stride == 32);
		I(rich_offset == 6 * 18 + 2 * 32);
		I(_rich_sector_stride == 32);

		if(rich) {
			I(reduced_id >= 0 &&
			  reduced_id < _rich_sector_stride);

			_id = rich_offset +
				continuous_sector * _rich_sector_stride +
				reduced_id;

			I(_id >= rich_offset && _id < nsupermodule);
			I(this->rich());
			I(this->sector() == continuous_sector);
			I(this->reduced_id() == reduced_id);
		}
		else {
			if(continuous_sector <
			   static_cast<int>(_nsector_lead_scintillator)) {
#if 0
				I(reduced_id >= 0 &&
				  reduced_id < _lead_scintillator_sector_stride);
#endif

				_id = lead_scintillator_offset +
					continuous_sector *
					_lead_scintillator_sector_stride +
					reduced_id;

#if 0
				I(_id >= 0 && _id < nsupermodule);
				I(this->lead_scintillator());
				I(this->sector() == continuous_sector);
				I(this->reduced_id() == reduced_id);
#endif
			}
			else {
				I(reduced_id >= 0 &&
				  reduced_id < _lead_glass_sector_stride);

				_id = lead_glass_offset +
					(continuous_sector - 6) *
					_lead_glass_sector_stride +
					reduced_id;

				I(_id >= 0 && _id < nsupermodule);
				I(this->lead_glass());
				I(this->sector() == continuous_sector);
				I(this->reduced_id() == reduced_id);
			}
		}
	}

	supermodule_id_t::
	supermodule_id_t(const int arm, const int sector,
					 const int reduced_id, const bool rich)
	{
		set_supermodule_id(arm, sector, reduced_id, rich);
	}

	supermodule_id_t::
	supermodule_id_t(const int continuous_sector,
					 const int reduced_id, const bool rich)
	{
		int arm;
		int sector;

		if(continuous_sector < 4) {
			arm = 0;
			sector = continuous_sector;
		}
		else {
			arm = 1;
			sector = continuous_sector < 6 ?
				continuous_sector - 2 : continuous_sector - 6;
		}
		set_supermodule_id(arm, sector, reduced_id, rich);
	}

	int supermodule_id_t::sector(void) const
	{
		if(lead_scintillator()) {
			const int ret = (_id - lead_scintillator_offset) /
				_lead_scintillator_sector_stride;

			I(ret >= 0 && ret < 6);

			return ret;
		}
		if(lead_glass()) {
			const int ret = 6 + (_id - lead_glass_offset) /
				_lead_glass_sector_stride;

			I(ret >= 6 && ret < 8);

			return ret;
		}

		return -1;
	}

	int supermodule_id_t::reduced_id(void) const
	{
		I(_lead_scintillator_sector_stride == 18);
		I(_lead_glass_sector_stride == 32);
		I(_rich_sector_stride == 32);

		if(lead_scintillator()) {
			const int ret = (_id - lead_scintillator_offset) %
				_lead_scintillator_sector_stride;

			I(ret >= 0 && ret < _lead_scintillator_sector_stride);

			return ret;
		}
		if(lead_glass()) {
			const int ret = (_id - lead_glass_offset) %
				_lead_glass_sector_stride;

			I(ret >= 0 && ret < _lead_glass_sector_stride);

			return ret;
		}

		return -1;
	}

	int supermodule_id_t::y(void) const
	{
		I(_lead_scintillator_sector_stride == 18);
		I(_lead_scintillator_z_stride == 6);
		I(_lead_glass_sector_stride == 32);
		I(_lead_glass_z_stride == 8);
		I(_rich_sector_stride == 32);
		I(_rich_z_stride == 8);

		if(lead_scintillator()) {
			const int ret = ((_id - lead_scintillator_offset) %
					_lead_scintillator_sector_stride) /
				_lead_scintillator_z_stride;

			I(ret >= 0 &&
			  ret < _lead_scintillator_sector_stride /
			  _lead_scintillator_z_stride);

			return ret;
		}
		if(lead_glass()) {
			const int ret = ((_id - lead_glass_offset) %
					_lead_glass_sector_stride) /
				_lead_glass_z_stride;

			I(ret >= 0 &&
			  ret < _lead_glass_sector_stride /
			  _lead_glass_z_stride);

			return ret;
		}

		return -1;
	}

	int supermodule_id_t::z(void) const
	{
		I(lead_scintillator_offset == 0);
		I(_lead_scintillator_z_stride == 6);
		I(lead_glass_offset == 6 * 18);
		I(_lead_glass_z_stride == 8);
		I(rich_offset == 6 * 18 + 2 * 32);
		I(_rich_z_stride == 8);

		if(lead_scintillator()) {
			const int ret = (_id - lead_scintillator_offset) %
				_lead_scintillator_z_stride;

			I(ret >= 0 && ret < _lead_scintillator_z_stride);

			return ret;
		}
		if(lead_glass()) {
			const int ret = (_id - lead_glass_offset) %
				_lead_glass_z_stride;

			I(ret >= 0 && ret < _lead_glass_z_stride);

			return ret;
		}

		return -1;
	}

	float supermodule_id_t::approximate_pseudorapidity(void) const
	{
		// Tracking and calorimetry pseudorapidity range
		static const float pseudorapidity_low = -0.35F;
		static const float pseudorapidity_high = 0.35F;

		if(rich()) {
			return NAN;
		}
		if(lead_scintillator()) {
			if(sector() < 4) {
				// Looking towards the interaction point, the PHENIX
				// coodinate system x-axis points towards east/out of
				// the east arm, y-axis points to the top, therefore
				// -z is to the right.
				return pseudorapidity_high -
					(z() + 0.5F) * (0.7F / 6.0F);
			}
			else {
				return pseudorapidity_low +
					(z() + 0.5F) * (0.7F / 6.0F);
			}
		}
		else {
			return pseudorapidity_low +
				(z() + 0.5F) * (0.7F / 8.0F);
		}
	}

	float supermodule_id_t::approximate_azimuth(void) const
	{
		// Calorimetry and infinite momentum tracking azimuth range
		static const float azimuth_west_low =
			-3.0F / 16.0F * (float)M_PI;
		static const float azimuth_east_high =
			-13.0F / 16.0F * (float)M_PI;

		if(rich()) {
			return NAN;
		}
		if(lead_scintillator()) {
			if(sector() < 4) {
				// Sectors 0, 1, 2, 3 = W0, W1, W2, W3
				return jet::angular_range_reduce(
					azimuth_west_low +
					sector() * (M_PI / 8.0) +
					(y() + 0.5F) * (M_PI / 24.0));
			}
			else {
				// Sectors 4, 5 = E2, E3
				return jet::angular_range_reduce(
					azimuth_east_high -
					((sector() - 2) * (M_PI / 8.0) +
					 (y() + 0.5F) * (M_PI / 24.0)));
			}
		}
		else {
			// Sectors 6, 7 = E0, E1
			return jet::angular_range_reduce(
				azimuth_east_high -
				((sector() - 6) * (M_PI / 8.0) +
				 (y() + 0.5F) * (M_PI / 32.0)));
		}
	}

#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 161)
#endif // __INTEL_COMPILER

	supermodule_id_t tower_id_t::supermodule(void) const
	{
		int id = -1;

		if(lead_scintillator()) {
			id = (y() / 12) * (_lead_scintillator_z_stride / 12) +
				(z() / 12);
		}
		if(lead_glass()) {
			id = (y() / 12) * (_lead_glass_z_stride / 12) +
				(z() / 12);
		}
		return supermodule_id_t(sector(), id, false);
	}

	std::string ert_hit_t::trigger_mode_str(void) const
	{
		switch(_trigger_mode) {
		case TRIGGER_MODE_4X4A:
			return "TRIGGER_MODE_4X4A";
			break;
		case TRIGGER_MODE_4X4B:
			return "TRIGGER_MODE_4X4B";
			break;
		case TRIGGER_MODE_4X4C:
			return "TRIGGER_MODE_4X4C";
			break;
		case TRIGGER_MODE_2X2:
			return "TRIGGER_MODE_2X2";
			break;
		case TRIGGER_MODE_RICH:
			return "TRIGGER_MODE_RICH";
			break;
		default:
			return "unknown trigger mode";
			break;
		}
	}

	void event_t::print(void) const
	{
		for(unsigned long i = 0; i < _track.size(); i++) {
#if 0
			if(_track[i].history().first_initial_state(_track) ==
			   _track.end())
#endif
				fprintf(stderr,
						"I*%c %3lu %3d %3d %8.3f %6.3f %6.3f %s\n",
						_track[i].history().is_initial_state() ? '*' : ' ',
						i,
						_track[i].history().parent_line_number().first,
						_track[i].history().parent_line_number().second,
						_track[i].momentum().perp(),
						_track[i].momentum().pseudorapidity(),
						_track[i].momentum().azimuth(),
						_track[i].symbol().c_str());
		}
	}

	event_t event_t::filter(bool (*track_filter)(const track_t),
							bool (*cluster_filter)(const cluster_t))
		const
	{
		event_t ret(_event_number, _raw_time);

#ifdef _OPENMP
#pragma omp sections
		{
#pragma omp section
#endif // _OPENMP
			for(unsigned long i = 0; i < _track.size(); i++)
				if(track_filter(_track[i]))
					ret._track.push_back(_track[i]);
#ifdef _OPENMP
#pragma omp section
#endif // _OPENMP
			for(unsigned long i = 0; i < _cluster.size(); i++)
				if(cluster_filter(_cluster[i]))
					ret._cluster.push_back(_cluster[i]);
#ifdef _OPENMP
#pragma omp section
#endif // _OPENMP
			ret._parton = _parton;
#ifdef _OPENMP
		}	// #pragma omp section
#endif // _OPENMP

		return ret;
	}

	event_t event_t::final_state_filter(void) const
	{
		event_t ret(_event_number, _raw_time);

#ifdef _OPENMP
#pragma omp sections
		{
#pragma omp section
#endif // _OPENMP
			for(unsigned long i = 0; i < _track.size(); i++)
				if(_track[i].is_final_state())
					ret._track.push_back(_track[i]);
#ifdef _OPENMP
#pragma omp section
#endif // _OPENMP
			for(unsigned long i = 0; i < _cluster.size(); i++)
				if(_cluster[i].is_final_state())
					ret._cluster.push_back(_cluster[i]);
#ifdef _OPENMP
#pragma omp section
#endif // _OPENMP
			ret._parton = _parton;
#ifdef _OPENMP
		}	// #pragma omp section
#endif // _OPENMP

		return ret;
	}

	event_t event_t::
	rectangular_filter(const std::pair<float, float>
					   pseudorapidity_range,
					   const std::pair<float, float>
					   azimuth_range) const
	{
		event_t ret(_event_number, _raw_time);

		if(azimuth_range.second >= azimuth_range.first) {
			// Range is azimuth_range.first <= azimuth <=
			// azimuth_range.second
			for(int i = 0; i < (int)_track.size(); i++) {
				const float pseudorapidity =
					_track[i].momentum().pseudorapidity();

				if(pseudorapidity >=
				   pseudorapidity_range.first &&
				   pseudorapidity <=
				   pseudorapidity_range.second) {
					const float azimuth =
						_track[i].momentum().azimuth();

					if(azimuth >= azimuth_range.first &&
					   azimuth <= azimuth_range.second)
						ret._track.push_back(_track[i]);
				}
			}
			for(int i = 0; i < (int)_cluster.size(); i++) {
				const float pseudorapidity =
					_cluster[i].momentum().pseudorapidity();

				if(pseudorapidity >=
				   pseudorapidity_range.first &&
				   pseudorapidity <=
				   pseudorapidity_range.second) {
					const float azimuth =
						_cluster[i].momentum().azimuth();

					if(azimuth >= azimuth_range.first &&
					   azimuth <= azimuth_range.second)
						ret._cluster.push_back(_cluster[i]);
				}
			}
		}
		else {
			// Range is azimuth_range.first <= azimuth <= M_PI or
			// -M_PI <= azimuth <= azimuth_range.second
			for(int i = 0; i < (int)_track.size(); i++) {
				const float pseudorapidity =
					_track[i].momentum().pseudorapidity();

				if(pseudorapidity >=
				   pseudorapidity_range.first &&
				   pseudorapidity <=
				   pseudorapidity_range.second) {
					const float azimuth =
						_track[i].momentum().azimuth();

					if(azimuth >= azimuth_range.first ||
					   azimuth <= azimuth_range.second)
						ret._track.push_back(_track[i]);
				}
			}
			for(int i = 0; i < (int)_cluster.size(); i++) {
				const float pseudorapidity =
					_cluster[i].momentum().pseudorapidity();

				if(pseudorapidity >=
				   pseudorapidity_range.first &&
				   pseudorapidity <=
				   pseudorapidity_range.second) {
					const float azimuth =
						_cluster[i].momentum().azimuth();

					if(azimuth >= azimuth_range.first ||
					   azimuth <= azimuth_range.second)
						ret._cluster.push_back(_cluster[i]);
				}
			}
		}
		ret._parton = _parton;

		return ret;
	}

	event_t event_t::tower_scale(const float tower_scale_scale[])
		const
	{
		event_t ret = *this;

		ret.cluster().clear();

		const unsigned long cluster_size = _cluster.size();

		ret._cluster.reserve(cluster_size);
		for(unsigned int i = 0; i < cluster_size; i++) {
			cluster_t cluster = _cluster[i];
			const int tower_id = cluster.central_tower_id();

			cluster.momentum().time() =
				tower_scale_scale[tower_id] *
				_cluster[i].momentum().time();
			ret._cluster.push_back(cluster);
		}
		ret._parton = _parton;

		return ret;
	}

	event_t event_t::tower_map_filter(const bool tower_map[]) const
	{
		event_t ret = *this;

		ret.cluster().clear();

		const unsigned long cluster_size = _cluster.size();

		ret._cluster.reserve(cluster_size);
		for(unsigned int i = 0; i < cluster_size; i++) {
			if(_cluster[i].tower_id_match(tower_map)) {
				I(tower_map[_cluster[i].central_tower_id()]);

				ret._cluster.push_back(_cluster[i]);
			}
		}
		ret._parton = _parton;

		return ret;
	}

#if 1
	void event_t::tower_diagnostics(const std::pair<float, float>
									pseudorapidity_range,
									const std::pair<float, float>
									azimuth_range, const float threshold) const
	{
		static std::vector<uint64_t> count(tower_id_t::ntower, 0);
		static std::vector<double> energy_count(tower_id_t::ntower, 0);
		static uint64_t call_count(0);

		if(azimuth_range.second >= azimuth_range.first) {
			// Range is azimuth_range.first <= azimuth <=
			// azimuth_range.second
			for(int i = 0; i < (int)_cluster.size(); i++) {
				const float pseudorapidity =
					_cluster[i].momentum().pseudorapidity();

				if(pseudorapidity >=
				   pseudorapidity_range.first &&
				   pseudorapidity <=
				   pseudorapidity_range.second) {
					const float azimuth =
						_cluster[i].momentum().azimuth();

					if(azimuth >= azimuth_range.first &&
					   azimuth <= azimuth_range.second) {
						for(int j = 0;
							j < (int)_cluster[i].tower_id().size();
							j++) {
							const int id = _cluster[i].tower_id()[j];

							I(id >= 0 && id < tower_id_t::ntower);

							count[id]++;
							energy_count[id] += _cluster[i].energy();
						}
					}
				}
			}
		}
		else {
			// Range is azimuth_range.first <= azimuth <= M_PI or
			// -M_PI <= azimuth <= azimuth_range.second
			for(int i = 0; i < (int)_cluster.size(); i++) {
				const float pseudorapidity =
					_cluster[i].momentum().pseudorapidity();

				if(pseudorapidity >=
				   pseudorapidity_range.first &&
				   pseudorapidity <=
				   pseudorapidity_range.second) {
					const float azimuth =
						_cluster[i].momentum().azimuth();

					if(azimuth >= azimuth_range.first ||
					   azimuth <= azimuth_range.second) {
						for(int j = 0;
							j < (int)_cluster[i].tower_id().size();
							j++) {
							const int id = _cluster[i].tower_id()[j];

							I(id >= 0 && id < tower_id_t::ntower);

							count[id]++;
							energy_count[id] += _cluster[i].energy();
						}
					}
				}
			}
		}

		if(call_count % 100 == 0) {
			float moment_1 = 0;
			float moment_2 = 0;
			int nonzero_count = 0;

			for(unsigned int i = 0; i < tower_id_t::ntower; i++) {
				if(count[i] > 0) {
					uint64_t c = count[i];

					moment_1 += c;
					moment_2 += c * c;
					nonzero_count++;
				}
			}
			moment_1 /= nonzero_count;
			moment_2 /= nonzero_count;

#if 0
			const float standard_deviation =
				sqrtf(moment_2 - moment_1 * moment_1);
#endif

			for(unsigned int i = 0; i < tower_id_t::ntower; i++) {
#if 0
				if(count[i] >
				   moment_1 + threshold * standard_deviation) {
#else
				if(std::isfinite(threshold) && count[i] > 0) {
#endif
					uint64_t c = count[i];
					double e = energy_count[i];

					std::cerr << i << '\t' << c << '\t' << e << '\t'
							  << e / c << std::endl;
				}
			}
			std::cerr << "------------------------------------------"
				"---------------------------" << std::endl;
		}
		call_count++;
	}
#endif

	float event_t::total_perp(void) const
	{
		const int track_size = _track.size();
		const int cluster_size = _cluster.size();
		float sum_track = 0;
		float sum_cluster = 0;

#ifdef _OPENMP
#pragma omp sections
		{
#pragma omp section
#pragma omp parallel for reduction(+: sum_track)
#endif // _OPENMP
			for(int i = 0; i < track_size; i++) {
				float perp = _track[i].momentum().perp();

				sum_track += std::isfinite(perp) ? perp : 0;
			}

#ifdef _OPENMP
#pragma omp section
#pragma omp parallel for reduction(+: sum_cluster)
#endif // _OPENMP
			for(int i = 0; i < cluster_size; i++) {
				float perp = _cluster[i].momentum().perp();

				sum_cluster += std::isfinite(perp) ? perp : 0;
			}
#ifdef _OPENMP
		}	// #pragma omp section
#endif // _OPENMP

		return sum_track + sum_cluster;
	}

	event_t event_t::operator+(const event_t &e) const
	{
		event_t ret = *this;

#ifdef _OPENMP
#pragma omp sections
		{
#pragma omp section
#endif // _OPENMP
			if(!(ret._centrality >= 0 &&
				 ret._centrality <= 100)) {
				// *this is a p + p event, e is presumed to be a heavy
				// ion event
				ret._vertex_height = e._vertex_height;
				ret._centrality = e._centrality;
				ret._reaction_plane = e._reaction_plane;
			}
			else {
				// *this is a heavy ion event, e is presumed to be a p
				// + p event
				ret._event_number = e._event_number;
				ret._raw_time = e._raw_time;
				ret._scale_down = e._scale_down;
				ret._ert_hit = e._ert_hit;
			}
#ifdef _OPENMP
#pragma omp section
#endif // _OPENMP
			for(int i = 0; i < static_cast<int>(e._track.size());
				i++) {
				ret._track.push_back(e._track[i]);
			}
#ifdef _OPENMP
#pragma omp section
#endif // _OPENMP
			for(int i = 0; i < static_cast<int>(e._cluster.size());
				i++) {
				ret._cluster.push_back(e._cluster[i]);
			}
#ifdef _OPENMP
#pragma omp section
#endif // _OPENMP
			for(int i = 0; i < static_cast<int>(e._parton.size());
				i++) {
				ret._parton.push_back(e._parton[i]);
			}
#ifdef _OPENMP
		}	// #pragma omp section
#endif // _OPENMP

		return ret;
	}

	event_t event_t::operator+=(const event_t &e)
	{
		*this = *this + e;

		return *this;
	}

	event_t event_t::vertex_align(const event_t &e) const
	{
		const float shift = static_cast<float>(e._vertex_height) -
			static_cast<float>(_vertex_height);
		event_t ret = *this;

		ret._track.clear();
		for(std::vector<track_t>::const_iterator iterator =
				_track.begin();
			iterator != _track.end(); iterator++) {
			track_t shifted_track = *iterator;

			shifted_track.pc1_crossing_height() += shift;
			ret._track.push_back(shifted_track);
		}
		ret._vertex_height += shift;

		return ret;
	}

	event_t event_t::rotate(const float azimuth) const
	{
		event_t ret = *this;

		ret._track.clear();
		for(std::vector<track_t>::const_iterator iterator =
				_track.begin();
			iterator != _track.end(); iterator++) {
			track_t rotated_track = *iterator;

			rotated_track.momentum().azimuth() =
				angular_range_reduce(
					iterator->momentum().azimuth() + azimuth);
			ret._track.push_back(rotated_track);
		}
		ret._cluster.clear();
		for(std::vector<cluster_t>::const_iterator iterator =
				_cluster.begin();
			iterator != _cluster.end(); iterator++) {
			cluster_t rotated_cluster = *iterator;

			rotated_cluster.momentum().azimuth() =
				angular_range_reduce(
					iterator->momentum().azimuth() + azimuth);
			ret._cluster.push_back(rotated_cluster);
		}

		return ret;
	}

#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif // __INTEL_COMPILER

}
