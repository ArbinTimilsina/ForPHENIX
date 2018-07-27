#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H
#include <jetevent/event.h>

namespace jet {

	const unsigned int particle_t::nuser_data = 16U;

    std::string particle_t::pdg_symbol(const int pdg_code)
    {
#ifndef PDG_NO_CHARGE_3
#define PDG_NO_CHARGE_3
#endif // PDG_NO_CHARGE_3
#include <table/pdg.h>
#undef PDG_NO_CHARGE_3

		const int *p = std::lower_bound(
			__pdg_code, __pdg_code + __npdg, pdg_code);
		if(p != __pdg_code + __npdg && *p == pdg_code) {
			const unsigned long index = p - __pdg_code;

			I(index < __npdg);

			return std::string(__pdg_symbol[index]);
		}
		else if(pdg_code == PDG_CODE_UNDEFINED) {
			return "(undefined)";
		}
		else if(pdg_code == PDG_CODE_INDETERMINATE) {
			return "(indeterminate)";
		}
		else if(pdg_code >= PDG_CODE_JET) {
			char buf[32];

			snprintf(buf, 32, "(jet_%d)", pdg_code - PDG_CODE_JET);

			return buf;
		}
		else {
			return "(invalid)";
		}
    }

    int particle_t::pdg_charge3(const int pdg_code)
    {
#ifndef PDG_NO_SYMBOL
#define PDG_NO_SYMBOL
#endif // PDG_NO_SYMBOL
#include <table/pdg.h>
#undef PDG_NO_SYMBOL

		const int *p = std::lower_bound(
			__pdg_code, __pdg_code + __npdg, pdg_code);
		if(p != __pdg_code + __npdg && *p == pdg_code) {
			const unsigned long index = p - __pdg_code;

			I(index < __npdg);

			return __pdg_charge_3[index];
		}
		else
			return 0;
    }

	bool particle_t::pdg_is_gluon(const int pdg_code)
	{
		return pdg_code == 21 || pdg_code == 9;
	}

	bool particle_t::pdg_is_photon(const int pdg_code)
	{
		return pdg_code == 22;
	}

	bool particle_t::pdg_is_z0(const int pdg_code)
	{
		return pdg_code == 23;
	}

	bool particle_t::pdg_is_quark(const int pdg_code)
	{
		const int pdg_code_particle_unexcited =
			std::abs(pdg_code) % 100000;

		return (pdg_code_particle_unexcited >= 1 &&
				pdg_code_particle_unexcited <= 8);
	}

	bool particle_t::pdg_is_diquark(const int pdg_code)
	{
		const int pdg_code_particle_unexcited =
			std::abs(pdg_code) % 100000;

		return (pdg_code_particle_unexcited >= 1103 &&
				pdg_code_particle_unexcited <= 8803 &&
				(pdg_code_particle_unexcited % 100) / 10 == 0);
	}

	bool particle_t::pdg_is_light_quark_or_diquark(const int pdg_code)
	{
		const int pdg_code_particle_unexcited =
			std::abs(pdg_code) % 100000;

		return ((pdg_code_particle_unexcited >= 1 &&
				 pdg_code_particle_unexcited <= 3) ||
				(pdg_code_particle_unexcited >= 1103 &&
				 pdg_code_particle_unexcited <= 8803 &&
				 (pdg_code_particle_unexcited % 100) / 10 == 0 &&
				 ((pdg_code_particle_unexcited / 1000 >= 1 &&
				   pdg_code_particle_unexcited / 1000 <= 3) ||
				  ((pdg_code_particle_unexcited % 1000) / 100 >= 1 &&
				   (pdg_code_particle_unexcited % 1000) / 100 <= 3))));
	}

	bool particle_t::pdg_is_charm_quark_or_diquark(const int pdg_code)
	{
		const int pdg_code_particle_unexcited =
			std::abs(pdg_code) % 100000;

		return (pdg_code_particle_unexcited == 4 ||
				(pdg_code_particle_unexcited >= 1103 &&
				 pdg_code_particle_unexcited <= 8803 &&
				 (pdg_code_particle_unexcited % 100) / 10 == 0 &&
				 (pdg_code_particle_unexcited / 1000 == 4 ||
				  (pdg_code_particle_unexcited % 1000) / 100 == 4)));
	}

	bool particle_t::pdg_is_bottom_quark_or_diquark(const int pdg_code)
	{
		const int pdg_code_particle_unexcited =
			std::abs(pdg_code) % 100000;

		return (pdg_code_particle_unexcited == 5 ||
				(pdg_code_particle_unexcited >= 1103 &&
				 pdg_code_particle_unexcited <= 8803 &&
				 (pdg_code_particle_unexcited % 100) / 10 == 0 &&
				 (pdg_code_particle_unexcited / 1000 == 5 ||
				  (pdg_code_particle_unexcited % 1000) / 100 == 5)));
	}

	bool particle_t::pdg_is_charm_meson(const int pdg_code)
	{
		const int pdg_code_particle_angular_momentum =
			std::abs(pdg_code) % 10000;

		return (pdg_code_particle_angular_momentum / 1000 == 0 &&
				(pdg_code_particle_angular_momentum % 1000) / 100 == 4 &&
				((pdg_code_particle_angular_momentum % 100) / 10 >= 1 &&
				 (pdg_code_particle_angular_momentum % 100) / 10 <= 3) &&
				(pdg_code_particle_angular_momentum % 10 == 1 ||
				 pdg_code_particle_angular_momentum % 10 == 3 ||
				 pdg_code_particle_angular_momentum % 10 == 5));
	}

	bool particle_t::pdg_is_bottom_meson(const int pdg_code)
	{
		const int pdg_code_particle_angular_momentum =
			std::abs(pdg_code) % 10000;

		return (pdg_code_particle_angular_momentum / 1000 == 0 &&
				(pdg_code_particle_angular_momentum % 1000) / 100 == 5 &&
				((pdg_code_particle_angular_momentum % 100) / 10 >= 1 &&
				 (pdg_code_particle_angular_momentum % 100) / 10 <= 4) &&
				(pdg_code_particle_angular_momentum % 10 == 1 ||
				 pdg_code_particle_angular_momentum % 10 == 3 ||
				 pdg_code_particle_angular_momentum % 10 == 5));
	}

	bool particle_t::pdg_is_neutral_pion(const int pdg_code)
	{
		const int pdg_code_particle_unexcited =
			std::abs(pdg_code) % 100000;

		return pdg_code_particle_unexcited == 111;
	}

    std::string particle_t::symbol(void) const
    {
		return pdg_symbol(_pdg_code);
    }

    int particle_t::charge3(void) const
    {
		return pdg_charge3(_pdg_code);
    }

#ifndef __CINT__
    rational_t particle_t::charge(void) const
    {
		const int q3 = charge3();

		return rational_t(q3, 3);
    }
#endif // __CINT__

	uint16_t particle_t::compact_pdg_code(void) const
	{
#ifndef PDG_NO_CHARGE_3
#define PDG_NO_CHARGE_3
#endif // PDG_NO_CHARGE_3
#ifndef PDG_NO_SYMBOL
#define PDG_NO_SYMBOL
#endif // PDG_NO_SYMBOL
#include <table/pdg.h>
#undef PDG_NO_SYMBOL
#undef PDG_NO_CHARGE_3

		const int *p = std::lower_bound(
			__pdg_code, __pdg_code + __npdg, _pdg_code);
		if(p != __pdg_code + __npdg && *p == _pdg_code) {
			const unsigned long index = p - __pdg_code;

			I(index < __npdg);

			return (uint16_t)index;
		}
		else if(_pdg_code >= PDG_CODE_JET &&
				_pdg_code < PDG_CODE_JET + 1920) {
			return (USHRT_MAX - 2048) + (_pdg_code - PDG_CODE_JET);
		}
		else {
			return USHRT_MAX;
		}
	}

	void particle_t::
	set_compact_pdg_code(const uint16_t compact_pdg_code)
	{
#ifndef PDG_NO_CHARGE_3
#define PDG_NO_CHARGE_3
#endif // PDG_NO_CHARGE_3
#ifndef PDG_NO_SYMBOL
#define PDG_NO_SYMBOL
#endif // PDG_NO_SYMBOL
#include <table/pdg.h>
#undef PDG_NO_SYMBOL
#undef PDG_NO_CHARGE_3

		if(compact_pdg_code < __npdg) {
			_pdg_code = __pdg_code[compact_pdg_code];
		}
		else if(compact_pdg_code >= (USHRT_MAX - 2048) &&
				compact_pdg_code < (USHRT_MAX - 2048) + 1920) {
			_pdg_code = PDG_CODE_JET +
				(compact_pdg_code - (USHRT_MAX - 2048));
		}
		else {
			_pdg_code = PDG_CODE_UNDEFINED;
		}
	}

    particle_t particle_t::snowmass_rotate(const float angle)
    {
		const float eta = _momentum[2];
		const float phi = _momentum[3];
#ifdef _GNU_SOURCE
		float sin_angle;
		float cos_angle;

		sincosf(angle, &sin_angle, &cos_angle);
#else // _GNU_SOURCE
		const float sin_angle = sinf(angle);
		const float cos_angle = cosf(angle);
#endif // _GNU_SOURCE

		_momentum[2] = cos_angle * eta - sin_angle * phi;
		_momentum[3] = sin_angle * eta + cos_angle * phi;

		return *this;
    }

    particle_t particle_t::operator+=(snowmass_vector_t v)
    {
		_momentum[2] += v[2];
		_momentum[3] = angular_range_reduce(_momentum[3] + v[3]);

		return *this;
    }

    particle_t particle_t::operator-=(snowmass_vector_t v)
    {
		_momentum[2] -= v[2];
		_momentum[3] = angular_range_reduce(_momentum[3] - v[3]);

		return *this;
    }

	bool particle_t::first_history_contain_bottom_quark_or_diquark(
		const std::vector<jet::track_t> particle_list) const
	{
		return first_history_contain(
			std::mem_fun_ref(&jet::particle_t::is_bottom_quark_or_diquark),
			particle_list);
	}

	bool particle_t::first_history_contain_charm_quark_or_diquark(
		const std::vector<jet::track_t> particle_list) const
	{
		return first_history_contain(
			std::mem_fun_ref(&jet::particle_t::is_charm_quark_or_diquark),
			particle_list);
	}

    /////////////////////////////////////////////////////////////////

    track_t::track_t(const particle_t particle, const int quality,
					 const int cluster_id,
					 const spherical_point_t dc_crossing,
					 const float pc1_crossing_height,
					 const rich_ring_t rich_ring,
					 const cylindrical_point_t
					 pc2_sigma_displacement,
					 const angular_slope_t tec_sigma_displacement,
					 const cylindrical_point_t
					 pc3_sigma_displacement,
					 const cylindrical_point_t
					 cluster_sigma_displacement)
		: particle_t(particle), _quality(quality),
		  _cluster_id(cluster_id), _dc_crossing(dc_crossing),
		  _pc1_crossing_height(pc1_crossing_height),
		  _rich_ring(rich_ring),
		  _pc2_sigma_displacement(pc2_sigma_displacement),
		  _tec_sigma_displacement(tec_sigma_displacement),
		  _pc3_sigma_displacement(pc3_sigma_displacement),
		  _cluster_sigma_displacement(cluster_sigma_displacement)
    {
		const int q3 = particle.charge3();

		_charge = q3 > 0 ? 1 : q3 < 0 ? -1 : 0;
    }
    track_t::track_t(const float energy, const float perp,
					 const float pseudorapidity, const float azimuth,
					 const int charge, const int quality,
					 const int cluster_id,
					 const float dc_crossing_polar_angle,
					 const float dc_crossing_azimuth,
					 const float pc1_crossing_height,
					 const int rich_ring_nnormal_area_phototube,
					 const float rich_ring_nnormal_area_photoelectron,
					 const int rich_ring_nwide_area_phototube,
					 const float rich_ring_nwide_area_photoelectron,
					 const float rich_ring_unitful_chi_square,
					 const float rich_ring_displacement,
					 const float pc2_sigma_displacement_azimuth,
					 const float pc2_sigma_displacement_height,
					 const float tec_sigma_displacement_azimuth,
					 const float tec_sigma_displacement_inclination,
					 const float pc3_sigma_displacement_azimuth,
					 const float pc3_sigma_displacement_height,
					 const float cluster_sigma_displacement_azimuth,
					 const float cluster_sigma_displacement_height)
		: particle_t(energy, perp, pseudorapidity, azimuth),
		  _charge(charge), _quality(quality),
		  _cluster_id(cluster_id),
		  _dc_crossing(dc_crossing_polar_angle, dc_crossing_azimuth),
		  _pc1_crossing_height(pc1_crossing_height),
		  _rich_ring(
				rich_ring_nnormal_area_phototube,
				rich_ring_nnormal_area_photoelectron,
				rich_ring_nwide_area_phototube,
				rich_ring_nwide_area_photoelectron,
				rich_ring_unitful_chi_square,
				rich_ring_displacement),
		  _pc2_sigma_displacement(
				pc2_sigma_displacement_azimuth,
				pc2_sigma_displacement_height),
		  _tec_sigma_displacement(
				tec_sigma_displacement_azimuth,
				tec_sigma_displacement_inclination),
		  _pc3_sigma_displacement(
				pc3_sigma_displacement_azimuth,
				pc3_sigma_displacement_height),
		  _cluster_sigma_displacement(
				cluster_sigma_displacement_azimuth,
				cluster_sigma_displacement_height)
    {
    }
    track_t::track_t(const float momentum[4], const int charge,
					 const int quality, const int cluster_id,
					 const float dc_crossing[2],
					 const float pc1_crossing_height,
					 const int rich_ring_count[4],
					 const float rich_ring_measure[2],
					 const float pc2_sigma_displacement[2],
					 const float tec_sigma_displacement[2],
					 const float pc3_sigma_displacement[2],
					 const float cluster_sigma_displacement[2])
		: particle_t(momentum), _charge(charge), _quality(quality),
		  _cluster_id(cluster_id), _dc_crossing(dc_crossing),
		  _pc1_crossing_height(pc1_crossing_height),
		  _rich_ring(rich_ring_count, rich_ring_measure),
		  _pc2_sigma_displacement(pc2_sigma_displacement),
		  _tec_sigma_displacement(tec_sigma_displacement),
		  _pc3_sigma_displacement(pc3_sigma_displacement),
		  _cluster_sigma_displacement(cluster_sigma_displacement)
    {
    }

	const unsigned int phenix_calorimetry_t::
	_lead_scintillator_sector_stride =
		_lead_scintillator_y_stride * _lead_scintillator_z_stride;

	const unsigned int phenix_calorimetry_t::
	_lead_glass_sector_stride =
		_lead_glass_y_stride * _lead_glass_z_stride;

	const unsigned int tower_id_t::lead_scintillator_offset = 0;

	const unsigned int tower_id_t::lead_glass_offset =
		_nsector_lead_scintillator *
		_lead_scintillator_sector_stride;

	const unsigned int tower_id_t::ntower =
		_nsector_lead_scintillator *
		_lead_scintillator_sector_stride +
		_nsector_lead_glass * _lead_glass_sector_stride;

	tower_id_t::tower_id_t(const int arm, const int sector,
						   const int y, const int z)
	{
		I(arm >= 0 && arm < 2);
		I(sector >= 0 && sector < 4);

		const int continuous_sector =
			arm == 0 ? sector :
			sector < 2 ? sector + 6 : sector + 2;

		I(continuous_sector >= 0 && continuous_sector < 8);

		if(continuous_sector < (int)_nsector_lead_scintillator) {
			I(z >= 0 && z < _lead_scintillator_z_stride);
			I(y >= 0 && y * _lead_scintillator_z_stride <
			  _lead_scintillator_sector_stride);

			_id = lead_scintillator_offset +
				continuous_sector *
				_lead_scintillator_sector_stride +
				y * _lead_scintillator_z_stride + z;

			I(_id >= 0 && _id < ntower);
			I(this->lead_scintillator());
			//I(this->sector() == continuous_sector);
			I(this->y() == y);
			I(this->z() == z);
		}
		else {
			I(z >= 0 && z < _lead_glass_z_stride);
			I(y >= 0 && y * _lead_glass_z_stride <
			  _lead_glass_sector_stride);

			_id = lead_glass_offset +
				(continuous_sector - 6) *
				_lead_glass_sector_stride +
				y * _lead_glass_z_stride + z;

			I(_id >= 0 && _id < ntower);
			I(this->lead_glass());
			//I(this->sector() == continuous_sector);
			I(this->y() == y);
			I(this->z() == z);
		}
	}

	int tower_id_t::sector(void) const
	{
		if(lead_scintillator()) {
			return (_id - lead_scintillator_offset) /
				_lead_scintillator_sector_stride;
		}
		if(lead_glass()) {
			return _nsector_lead_scintillator +
				(_id - lead_glass_offset) /
				_lead_glass_sector_stride;
		}

		return -1;
	}

	int tower_id_t::y(void) const
	{
		if(lead_scintillator()) {
			return ((_id - lead_scintillator_offset) %
					_lead_scintillator_sector_stride) /
				_lead_scintillator_z_stride;
		}
		if(lead_glass()) {
			return ((_id - lead_glass_offset) %
					_lead_glass_sector_stride) /
				_lead_glass_z_stride;
		}

		return -1;
	}

	int tower_id_t::z(void) const
	{
		if(lead_scintillator()) {
			return (_id - lead_scintillator_offset) %
				_lead_scintillator_z_stride;
		}
		if(lead_glass()) {
			return (_id - lead_glass_offset) %
				_lead_glass_z_stride;
		}

		return -1;
	}

    std::string nucleus_t::element_symbol(void) const
    {
		// Invalid A-Z range
		if(_nucleon_number < _proton_number)
			return "";

#include <table/iupac.h>

		I(__niupac >= 99);

		// Special non-IUPAC symbols for elementary particle and heavy
		// ion physics
		if(_proton_number == 0 && _nucleon_number == 1)
			return "n";					// Neutron
		else if(_proton_number == 1 && _nucleon_number <= 3) {
			switch(_nucleon_number) {
			case 1:	return "p"; break;	// Proton
			case 2:	return "d"; break;	// Deuterium
			case 3:	return "t"; break;	// Tritium
			}
		}

		if(_proton_number > 0 && _proton_number <= __niupac)
			// IUPAC element symbols
			return __iupac_symbol[_proton_number - 1];

		if(_proton_number > __niupac && _proton_number <= 999) {
			// IUPAC transuranic systematic symbols
			const int z0 = _proton_number / 100;
			const int z1 = (_proton_number / 10) % 10;
			const int z2 = _proton_number % 100;
			char buf[4];

			buf[0] = (char)toupper(__iupac_systematic_symbol_rule[z0]);
			buf[1] = __iupac_systematic_symbol_rule[z1];
			buf[2] = __iupac_systematic_symbol_rule[z2];
			buf[3] = '\0';

			return buf;
		}

		return "";
    }

    std::string nucleus_t::isotope_symbol(void) const
    {
		std::ostringstream stream;

		stream << element_symbol() << _nucleon_number;

		return stream.str();
    }

    nucleus_t proton = nucleus_t(1, 1);
    nucleus_t deuteron = nucleus_t(1, 2);
    nucleus_t copper_63 = nucleus_t(29, 63);
    nucleus_t gold_197 = nucleus_t(79, 197);
    nucleus_t lead_208 = nucleus_t(82, 208);
    nucleus_t uranium_238 = nucleus_t(92, 238);
	double rhic_sqrt_s_nn = 200;
	double lhc_sqrt_s_nn = 5520;
	double tevatron_sqrt_s = 1960;
	double lhc_sqrt_s = 14000;

    std::string collisional_system_t::name(void) const
    {
		std::ostringstream stream;
		std::string projectile_name = _projectile.name();
		std::string target_name = _target.name();
		std::string s_symbol = is_hadronic() ? "s" : "s_{NN}";

		stream << "\\mathrm{" << projectile_name
			   << "} + \\mathrm{" << target_name
			   << " \\:\\text{at}\\: \\sqrt{" << s_symbol
			   << "} = " << _cms_energy << "\\,\\mathrm{GeV}";

		return stream.str();
    }

}
