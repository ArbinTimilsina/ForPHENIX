// -*- mode: c++; -*-

#ifndef JETEVENT_PARTICLE_H_
#define JETEVENT_PARTICLE_H_

/////////////////////////////////////////////////////////////////////

// LinkDef files generated by rootcint do not include config.h
#if defined(__CINT__) || defined(G__DICTIONARY)
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H
#endif // defined(__CINT__) || defined(G__DICTIONARY)

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <jetbase/dbc.h>
#include <jetbase/num.h>
#include <jetbase/specfunc.h>
#include <jetbase/geometry.h>
#include <jetevent/snowmass.h>
#include <jetevent/lorentz.h>

namespace jet {

	class particle_t;
	class event_file_t;
	class reconstruction_filtering_t;

	class particle_history_t {
	private:
		std::pair<int, int> _parent_line_number;
		std::pair<int, int> _daughter_line_number;
	public:
		particle_history_t(void)
		{
		}
		particle_history_t(const int parent_line_number,
						   const int first_daughter_line_number,
						   const int second_daughter_line_number)
			: _parent_line_number(parent_line_number, -1),
			  _daughter_line_number(first_daughter_line_number,
									second_daughter_line_number)
		{
		}
		particle_history_t(const int first_parent_line_number,
						   const int second_parent_line_number,
						   const int first_daughter_line_number,
						   const int second_daughter_line_number)
			: _parent_line_number(first_parent_line_number,
								  second_parent_line_number),
			  _daughter_line_number(first_daughter_line_number,
									second_daughter_line_number)
		{
		}
		inline std::pair<int, int> parent_line_number(void) const
		{
			return _parent_line_number;
		}
		inline std::pair<int, int> &parent_line_number(void)
		{
			return _parent_line_number;
		}
		inline std::pair<int, int> daughter_line_number(void) const
		{
			return _daughter_line_number;
		}
		inline std::pair<int, int> &daughter_line_number(void)
		{
			return _daughter_line_number;
		}
		inline bool is_initial_state(void) const
		{
			return (_parent_line_number.first < 0 &&
					_parent_line_number.second < 0);
		}
		template<typename physical_particle_list_t>
		typename physical_particle_list_t::const_iterator
		first_initial_state(const physical_particle_list_t &
							particle_list) const
		{
			// The recursion accounts for three different cases:
			if(_parent_line_number.first < 0)
				// (a) the particle does not have any initial state
				// (return particle_list.end()),
				return particle_list.end();
			else if(particle_list[_parent_line_number.first].
					history().is_initial_state())
				// (b) the immediately previous particle is an initial
				// state (return it),
				return particle_list.begin() +
					_parent_line_number.first;
			else
				// (c) the particle is deeper in an ancestry tree
				// (continue by recursion on the tree).
				return particle_list[_parent_line_number.first].
					history().first_initial_state(particle_list);
		}
		// FIXME: full initial_state() has to be implemented
		friend class event_file_t;
	};

	/////////////////////////////////////////////////////////////////

	// (Monte-Carlo) particles, tracks and calorimeter clusters

	class track_t;

	// General particle class

#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 810)
#endif // __INTEL_COMPILER
	class particle_t {
	public:
		static const unsigned int nuser_data;
	protected:
		int _status_code;
		int _pdg_code;
		int64_t _user_data[16];
		snowmass_vector_t _momentum;
		lorentz_vector_t _vertex;
		particle_history_t _history;
		float _proper_lifetime;
		snowmass_vector_t _background_momentum;
	public:
		// GEANT particles plus unstable states useful for tagging
		enum {
			PDG_CODE_DOWN_QUARK = 1,
			PDG_CODE_DOWN_ANTIQUARK = -1,
			PDG_CODE_UP_QUARK = 2,
			PDG_CODE_UP_ANTIQUARK = -2,
			PDG_CODE_STRANGE_QUARK = 3,
			PDG_CODE_STRANGE_ANTIQUARK = -3,
			PDG_CODE_CHARM_QUARK = 4,
			PDG_CODE_CHARM_ANTIQUARK = -4,
			PDG_CODE_BOTTOM_QUARK = 5,
			PDG_CODE_BOTTOM_ANTIQUARK = -5,
			PDG_CODE_TOP_QUARK = 6,
			PDG_CODE_TOP_ANTIQUARK = -6,
			PDG_CODE_ELECTRON = 11,
			PDG_CODE_POSITRON = -11,
			PDG_CODE_ELECTRON_NEUTRINO = 12,
			PDG_CODE_ELECTRON_ANTINEUTRINO = -12,
			PDG_CODE_MUON = 13,
			PDG_CODE_ANTIMUON = -13,
			PDG_CODE_MUON_NEUTRINO = 14,
			PDG_CODE_MUON_ANTINEUTRINO = -14,
			PDG_CODE_TAU = 15,
			PDG_CODE_ANTITAU = -15,
			PDG_CODE_TAU_NEUTRINO = 16,
			PDG_CODE_TAU_ANTINEUTRINO = -16,
			PDG_CODE_GLUON = 21,
			PDG_CODE_PHOTON = 22,
			PDG_CODE_Z_0 = 23,
			PDG_CODE_W_PLUS = 24,
			PDG_CODE_W_MINUS = -24,
			PDG_CODE_GRAVITON = 39,
			PDG_CODE_PION_0 = 111,
			PDG_CODE_PION_PLUS = 211,
			PDG_CODE_PION_MINUS = -211,
			PDG_CODE_RHO_770_0 = 113,
			PDG_CODE_RHO_770_PLUS = 213,
			PDG_CODE_ANTIRHO_770_MINUS = -213,
			PDG_CODE_ETA = 221,
			PDG_CODE_OMEGA_782 = 223,
			PDG_CODE_KAON_0_LONG = 130,
			PDG_CODE_KAON_0_SHORT = 310,
			PDG_CODE_KAON_0 = 311,
			PDG_CODE_KAON_PLUS = 321,
			PDG_CODE_KAON_MINUS = -321,
			PDG_CODE_J_PSI = 443,
			PDG_CODE_PROTON = 2212,
			PDG_CODE_ANTIPROTON = -2212,
			PDG_CODE_NEUTRON = 2112,
			PDG_CODE_ANTINEUTRON = -2112,
			PDG_CODE_LAMBDA = 3122,
			PDG_CODE_ANTILAMBDA = -3122,
			PDG_CODE_SIGMA_PLUS = 3222,
			PDG_CODE_ANTISIGMA_MINUS = -3222,
			PDG_CODE_SIGMA_0 = 3212,
			PDG_CODE_ANTISIGMA_0 = -3212,
			PDG_CODE_SIGMA_MINUS = 3112,
			PDG_CODE_ANTISIGMA_PLUS = -3112,
			PDG_CODE_XI_0 = 3322,
			PDG_CODE_ANTIXI_0 = -3322,
			PDG_CODE_XI_MINUS = 3322,
			PDG_CODE_ANTIXI_PLUS = -3322,
			PDG_CODE_OMEGA_MINUS = 3334,
			PDG_CODE_ANTIOMEGA_PLUS = -3334
		};
		// Special purpose codes
		enum {
			PDG_CODE_UNDEFINED = INT_MIN,
			PDG_CODE_INDETERMINATE,
			PDG_CODE_INFRARED,
			PDG_CODE_JET = (1 << 24)
		};
		static std::string pdg_symbol(const int pdg_code);
		static int pdg_charge3(const int pdg_code);
		static bool pdg_is_gluon(const int pdg_code);
		static bool pdg_is_photon(const int pdg_code);
		static bool pdg_is_z0(const int pdg_code);
		static bool pdg_is_quark(const int pdg_code);
		static bool pdg_is_diquark(const int pdg_code);
		static bool pdg_is_light_quark_or_diquark(const int pdg_code);
		static bool pdg_is_charm_quark_or_diquark(const int pdg_code);
		static bool pdg_is_bottom_quark_or_diquark(const int pdg_code);
		static bool pdg_is_charm_meson(const int pdg_code);
		static bool pdg_is_bottom_meson(const int pdg_code);
		static bool pdg_is_neutral_pion(const int pdg_code);
		particle_t(const snowmass_vector_t momentum =
				   snowmass_vector_t(NAN, NAN, NAN, NAN),
				   const int pdg_code = PDG_CODE_UNDEFINED,
				   const int status_code = -1,
				   const lorentz_vector_t vertex =
				   lorentz_vector_t(NAN, NAN, NAN, NAN),
				   const particle_history_t history =
				   particle_history_t(-1, -1, -1, -1),
				   const float proper_lifetime = NAN,
				   const snowmass_vector_t background_momentum =
				   snowmass_vector_t(NAN, NAN, NAN, NAN))
			: _status_code(status_code), _pdg_code(pdg_code),
			  _momentum(momentum),
			  _vertex(vertex), _history(history),
			  _proper_lifetime(proper_lifetime),
			  _background_momentum(background_momentum)
		{
		}
		particle_t(const float energy, const float perp,
				   const float pseudorapidity, const float azimuth,
				   const int pdg_code = PDG_CODE_UNDEFINED,
				   const int status_code = -1,
				   const lorentz_vector_t vertex =
				   lorentz_vector_t(NAN, NAN, NAN, NAN),
				   const particle_history_t history =
				   particle_history_t(-1, -1, -1, -1),
				   const float proper_lifetime = NAN,
				   const snowmass_vector_t background_momentum =
				   snowmass_vector_t(NAN, NAN, NAN, NAN))
			: _status_code(status_code), _pdg_code(pdg_code),
			  _momentum(energy, perp, pseudorapidity, azimuth),
			  _vertex(vertex), _history(history),
			  _proper_lifetime(proper_lifetime),
			  _background_momentum(background_momentum)
		{
		}
		particle_t(const float momentum[4],
				   const int pdg_code = PDG_CODE_UNDEFINED,
				   const int status_code = 1,
				   const lorentz_vector_t vertex =
				   lorentz_vector_t(NAN, NAN, NAN, NAN),
				   const particle_history_t history =
				   particle_history_t(-1, -1, -1, -1),
				   const float proper_lifetime = NAN,
				   const snowmass_vector_t background_momentum =
				   snowmass_vector_t(NAN, NAN, NAN, NAN))
			: _status_code(status_code), _pdg_code(pdg_code),
			  _momentum(momentum),
			  _vertex(vertex), _history(history),
			  _proper_lifetime(proper_lifetime),
			  _background_momentum(background_momentum)
		{
		}
		inline snowmass_vector_t &momentum(void)
		{
			return _momentum;
		}
		inline snowmass_vector_t momentum(void) const
		{
			return _momentum;
		}
		inline float &energy(void)
		{
			return _momentum.time();
		}
		inline float energy(void) const
		{
			return _momentum.time();
		}
		inline float mass(void) const
		{
			return _momentum.mass();
		}
		inline float &pseudorapidity(void)
		{
			return _momentum.pseudorapidity();
		}
		inline float pseudorapidity(void) const
		{
			return _momentum.pseudorapidity();
		}
		inline float &azimuth(void)
		{
			return _momentum.azimuth();
		}
		inline float azimuth(void) const
		{
			return _momentum.azimuth();
		}
		inline int &status_code(void)
		{
			return _status_code;
		}
		inline int status_code(void) const
		{
			return _status_code;
		}
		inline int &pdg_code(void)
		{
			return _pdg_code;
		}
		inline int pdg_code(void) const
		{
			return _pdg_code;
		}
		// Explicit const for std::const_mem_fun_ref_t
		inline int get_pdg_code(void) const
		{
			return _pdg_code;
		}
		inline int64_t *user_data(void)
		{
			return _user_data;
		}
		inline const int64_t *user_data(void) const
		{
			return _user_data;
		}
		inline lorentz_vector_t vertex(void) const
		{
			return _vertex;
		}
		inline lorentz_vector_t &vertex(void)
		{
			return _vertex;
		}
		inline particle_history_t history(void) const
		{
			return _history;
		}
		inline particle_history_t &history(void)
		{
			return _history;
		}
		inline float proper_lifetime(void) const
		{
			return _proper_lifetime;
		}
		inline float &proper_lifetime(void)
		{
			return _proper_lifetime;
		}
		inline snowmass_vector_t background_momentum(void) const
		{
			return _background_momentum;
		}
		inline snowmass_vector_t &background_momentum(void)
		{
			return _background_momentum;
		}
		inline bool is_final_state(void) const
		{
			// This assumes Pythia and HERWIG convention
			return _status_code == 1;
		}
		std::string symbol(void) const;
		int charge3(void) const;
#ifndef __CINT__
		rational_t charge(void) const;
#endif // __CINT__
		inline bool is_gluon(void) const
		{
			return pdg_is_gluon(_pdg_code);
		}
		inline bool is_photon(void) const
		{
			return pdg_is_photon(_pdg_code);
		}
		inline bool is_z0(void) const
		{
			return pdg_is_z0(_pdg_code);
		}
		inline bool is_quark(void) const
		{
			return pdg_is_quark(_pdg_code);
		}
		inline bool is_diquark(void) const
		{
			return pdg_is_diquark(_pdg_code);
		}
		inline bool is_neutral_pion(void) const
		{
			return pdg_is_neutral_pion(_pdg_code);
		}
		inline bool is_light_quark_or_diquark(void) const
		{
			return pdg_is_light_quark_or_diquark(_pdg_code);
		}
		inline bool is_charm_quark_or_diquark(void) const
		{
			return pdg_is_charm_quark_or_diquark(_pdg_code);
		}
		inline bool is_bottom_quark_or_diquark(void) const
		{
			return pdg_is_bottom_quark_or_diquark(_pdg_code);
		}
		inline bool is_charm_meson(void) const
		{
			return pdg_is_charm_meson(_pdg_code);
		}
		inline bool is_bottom_meson(void) const
		{
			return pdg_is_bottom_meson(_pdg_code);
		}
		// Compact storage
		inline uint8_t compact_status_code(void) const
		{
			return _status_code < 0 ? 255U :
				static_cast<uint8_t>(_status_code);
		}
		inline void
		set_compact_status_code(const uint8_t compact_status_code)
		{
			_status_code = compact_status_code == 255U ? -1 :
				compact_status_code;
		}
		uint16_t compact_pdg_code(void) const;
		void set_compact_pdg_code(const uint16_t compact_pdg_code);
		particle_t snowmass_rotate(const float angle);
		particle_t operator+=(snowmass_vector_t v);
		particle_t operator-=(snowmass_vector_t v);
		inline bool operator<(const particle_t v) const
		{
			return _momentum < v._momentum;
		}
		template<typename physical_particle_list_t>
		inline bool first_history_contain(
			std::const_mem_fun_ref_t<bool, particle_t> test,
			const physical_particle_list_t particle_list) const
		{
			if(_history.parent_line_number().first < 0) {
				return test(*this);
			}
			else {
				return test(*this) ||
					particle_list[_history.parent_line_number().first].
					first_history_contain(test, particle_list);
			}
		}
#ifndef __CINT__
		bool first_history_contain_bottom_quark_or_diquark(
			const std::vector<jet::track_t> particle_list) const;
		bool first_history_contain_charm_quark_or_diquark(
			const std::vector<jet::track_t> particle_list) const;
#endif // __CINT__
		virtual ~particle_t(void)
		{
		}
		friend class event_file_t;
	};
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif // __INTEL_COMPILER

	// PHENIX central arms track

	class rich_ring_count_t {
	private:
		int _nphototube;
		float _nphotoelectron;
	public:
		rich_ring_count_t(void)
		{
		}
		rich_ring_count_t(const int nphototube,
						  const float nphotoelectron)
			: _nphototube(nphototube),
			  _nphotoelectron(nphotoelectron)
		{
		}
		int &nphototube(void)
		{
			return _nphototube;
		}
		int nphototube(void) const
		{
			return _nphototube;
		}
		float &nphotoelectron(void)
		{
			return _nphotoelectron;
		}
		float nphotoelectron(void) const
		{
			return _nphotoelectron;
		}
		friend class event_file_t;
	};

	class rich_ring_t {
	private:
		rich_ring_count_t _normal_area_count;
		rich_ring_count_t _wide_area_count;
		float _unitful_chi_square;
		float _sigma_displacement;
	public:
		inline rich_ring_t(const int nnormal_area_phototube = -1,
						   const float nnormal_area_photoelectron = NAN,
						   const int nwide_area_phototube = -1,
						   const float nwide_area_photoelectron = NAN,
						   const float unitful_chi_square = NAN,
						   const float sigma_displacement = NAN)
			: _normal_area_count(nnormal_area_phototube,
								 nnormal_area_photoelectron),
			  _wide_area_count(nwide_area_phototube,
							   nwide_area_photoelectron),
			  _unitful_chi_square(unitful_chi_square),
			  _sigma_displacement(sigma_displacement)
		{
		}
		inline rich_ring_t(const int count[2], const float measure[4])
			: _normal_area_count(count[0], measure[0]),
			  _wide_area_count(count[1], measure[1]),
			  _unitful_chi_square(measure[2]),
			  _sigma_displacement(measure[3])
		{
		}
		inline rich_ring_count_t &normal_area_count(void)
		{
			return _normal_area_count;
		}
		inline rich_ring_count_t normal_area_count(void) const
		{
			return _normal_area_count;
		}
		inline rich_ring_count_t &wide_area_count(void)
		{
			return _wide_area_count;
		}
		inline rich_ring_count_t wide_area_count(void) const
		{
			return _wide_area_count;
		}
		inline float &unitful_chi_square(void)
		{
			return _unitful_chi_square;
		}
		inline float unitful_chi_square(void) const
		{
			return _unitful_chi_square;
		}
		inline float &sigma_displacement(void)
		{
			return _sigma_displacement;
		}
		inline float sigma_displacement(void) const
		{
			return _sigma_displacement;
		}
		friend class event_file_t;
	};

	class track_t : public particle_t
	{
	public:
		enum {
			QUALITY_X1_USED = (1U << 0),
			QUALITY_X2_USED = (1U << 1),
			QUALITY_UV_FOUND = (1U << 2),
			QUALITY_UV_UNIQUE = (1U << 3),
			QUALITY_PC1_FOUND = (1U << 4),
			QUALITY_PC1_UNIQUE = (1U << 5)
		};
	private:
		int _charge;
		int _quality;
		int _cluster_id;
		// Matching variables, from inner to outer detectors
		spherical_point_t _dc_crossing;
		float _pc1_crossing_height;
		rich_ring_t _rich_ring;
		cylindrical_point_t _pc2_sigma_displacement;
		angular_slope_t _tec_sigma_displacement;
		cylindrical_point_t _pc3_sigma_displacement;
		cylindrical_point_t _cluster_sigma_displacement;
		std::vector<std::vector<jet::track_t>::const_iterator>
		_constituent;
	public:
		inline track_t(void)
		{
		}
		track_t(const particle_t particle,
				const int quality = 0,
				const int cluster_id = -1,
				const spherical_point_t dc_crossing =
				spherical_point_t(NAN, NAN),
				const float pc1_crossing_height = NAN,
				const rich_ring_t rich_ring = rich_ring_t(),
				const cylindrical_point_t pc2_sigma_displacement =
				cylindrical_point_t(NAN, NAN),
				const angular_slope_t tec_sigma_displacement =
				angular_slope_t(NAN, NAN),
				const cylindrical_point_t pc3_sigma_displacement =
				cylindrical_point_t(NAN, NAN),
				const cylindrical_point_t
				cluster_sigma_displacement =
				cylindrical_point_t(NAN, NAN));
		track_t(const float energy, const float perp,
				const float pseudorapidity, const float azimuth,
				const int charge = INT_MIN,
				const int quality = 0,
				const int cluster_id = -1,
				const float dc_crossing_polar_angle = NAN,
				const float dc_crossing_azimuth = NAN,
				const float pc1_crossing_height = NAN,
				const int rich_ring_nnormal_area_phototube = -1,
				const float rich_ring_nnormal_area_photoelectron = NAN,
				const int rich_ring_nwide_area_phototube = -1,
				const float rich_ring_nwide_area_photoelectron = NAN,
				const float rich_ring_unitful_chi_square = NAN,
				const float rich_ring_displacement = NAN,
				const float pc2_sigma_displacement_azimuth = NAN,
				const float pc2_sigma_displacement_height = NAN,
				const float tec_sigma_displacement_azimuth = NAN,
				const float tec_sigma_displacement_inclination = NAN,
				const float pc3_sigma_displacement_azimuth = NAN,
				const float pc3_sigma_displacement_height = NAN,
				const float cluster_sigma_displacement_azimuth = NAN,
				const float cluster_sigma_displacement_height = NAN);
		track_t(const float momentum[4], const int charge,
				const int quality, const int cluster_id,
				const float dc_crossing[2],
				const float pc1_crossing_height,
				const int rich_ring_count[4],
				const float rich_ring_measure[2],
				const float pc2_sigma_displacement[2],
				const float tec_sigma_displacement[2],
				const float pc3_sigma_displacement[2],
				const float cluster_sigma_displacement[2]);
		inline int &quality(void)
		{
			return _quality;
		}
		inline int quality(void) const
		{
			return _quality;
		}
		inline int &charge(void)
		{
			return _charge;
		}
		inline int charge(void) const
		{
			return _charge;
		}
		inline int &cluster_id(void)
		{
			return _cluster_id;
		}
		inline int cluster_id(void) const
		{
			return _cluster_id;
		}
		inline spherical_point_t &dc_crossing(void)
		{
			return _dc_crossing;
		}
		inline spherical_point_t dc_crossing(void) const
		{
			return _dc_crossing;
		}
		inline float &pc1_crossing_height(void)
		{
			return _pc1_crossing_height;
		}
		inline float pc1_crossing_height(void) const
		{
			return _pc1_crossing_height;
		}
		inline rich_ring_t &rich_ring(void)
		{
			return _rich_ring;
		}
		inline rich_ring_t rich_ring(void) const
		{
			return _rich_ring;
		}
		inline cylindrical_point_t &pc2_sigma_displacement(void)
		{
			return _pc2_sigma_displacement;
		}
		inline cylindrical_point_t pc2_sigma_displacement(void) const
		{
			return _pc2_sigma_displacement;
		}
		inline angular_slope_t &tec_sigma_displacement(void)
		{
			return _tec_sigma_displacement;
		}
		inline angular_slope_t tec_sigma_displacement(void) const
		{
			return _tec_sigma_displacement;
		}
		inline cylindrical_point_t &pc3_sigma_displacement(void)
		{
			return _pc3_sigma_displacement;
		}
		inline cylindrical_point_t pc3_sigma_displacement(void) const
		{
			return _pc3_sigma_displacement;
		}
		inline cylindrical_point_t &cluster_sigma_displacement(void)
		{
			return _cluster_sigma_displacement;
		}
		inline cylindrical_point_t cluster_sigma_displacement(void) const
		{
			return _cluster_sigma_displacement;
		}
		inline std::vector<std::vector<jet::track_t>::const_iterator> &
		constituent(void)
		{
			return _constituent;
		}
		inline std::vector<std::vector<jet::track_t>::const_iterator>
		constituent(void) const
		{
			return _constituent;
		}
		inline bool operator<(const track_t t) const
		{
			return _momentum < t._momentum;
		}
		friend class event_file_t;
		friend class reconstruction_filtering_t;
	};

	class phenix_calorimetry_t {
	protected:
		static const unsigned int _nsector_lead_scintillator = 6;
		static const unsigned int _nsector_lead_glass = 2;
		static const unsigned int _lead_scintillator_z_stride = 72;
		static const unsigned int _lead_scintillator_y_stride = 36;
		static const unsigned int _lead_glass_z_stride = 96;
		static const unsigned int _lead_glass_y_stride = 48;
		static const unsigned int _lead_scintillator_sector_stride;
		static const unsigned int _lead_glass_sector_stride;
		inline virtual ~phenix_calorimetry_t(void)
		{
		}
	};

	class supermodule_id_t;

	class tower_id_t : public phenix_calorimetry_t {
		// From the nDST classes (with a "iz" typo):
		//
		// towerid = (isect >= 6 ? 15552 + 4608 * (isect - 6) + 96 *
		// (4 * ismy + iz) + 6*ismz + iz : 2592 * isect + 72 * (12 *
		// ismy + iy) + 12 * ismz + iz);
	private:
		int _id;
	public:
		static const unsigned int lead_scintillator_offset;
		static const unsigned int lead_glass_offset;
		static const unsigned int ntower;
		inline tower_id_t(const int id = -1)
			: _id(id)
		{
		}
		tower_id_t(const int arm, const int sector,
				   const int y, const int z);
		inline operator int(void) const
		{
			return _id;
		}
		inline bool lead_scintillator(void) const
		{
			return (_id >= 0 &&
					_id < static_cast<int>(lead_glass_offset));
		}
		inline bool lead_glass(void) const
		{
			return (_id >= static_cast<int>(lead_glass_offset) &&
					_id < static_cast<int>(ntower));
		}
		int sector(void) const;
		int y(void) const;
		int z(void) const;
		inline bool operator==(const tower_id_t &id) const
		{
			return _id == id._id;
		}
		inline bool operator!=(const tower_id_t &id) const
		{
			return _id != id._id;
		}
#ifndef __CINT__
		supermodule_id_t supermodule(void) const;
		bool in(const supermodule_id_t &id) const;
#endif // __CINT__
		friend class event_file_t;
	};

	class tof_timing_t {
	private:
		float _mean;
		float _dispersion;
		float _min;
		float _max;
	public:
		tof_timing_t(const float mean = NAN,
					 const float dispersion = NAN,
					 const float min = NAN, const float max = NAN)
			: _mean(mean), _dispersion(dispersion),
			  _min(min), _max(max)
		{
		}
		inline float mean(void) const
		{
			return _mean;
		}
		inline float &mean(void)
		{
			return _mean;
		}
		inline float dispersion(void) const
		{
			return _dispersion;
		}
		inline float &dispersion(void)
		{
			return _dispersion;
		}
		inline float min(void) const
		{
			return _min;
		}
		inline float &min(void)
		{
			return _min;
		}
		inline float max(void) const
		{
			return _max;
		}
		inline float &max(void)
		{
			return _max;
		}
	};

	// PHENIX EMC clusters
	class cluster_t : public particle_t {
	private:
		int _id;
		std::vector<jet::tower_id_t> _tower_id;
		std::vector<float> _partial_energy_sum;
		tower_id_t _central_tower_id;
		float _electromagnetic_chi_square;
		tof_timing_t _tof_timing;
		float _incident_angle;
		float _corrected_dispersion_y;
		float _corrected_dispersion_z;
	public:
		cluster_t(void)
		{
		}
		cluster_t(const particle_t particle, const int id,
				  const std::vector<tower_id_t> tower_id =
				  std::vector<tower_id_t>(),
				  const std::vector<float> partial_energy_sum =
				  std::vector<float>(),
				  const float electromagnetic_chi_square = NAN,
				  const tof_timing_t tof_timing =
				  tof_timing_t(NAN, NAN, NAN, NAN),
				  const float incident_angle = NAN,
				  const float corrected_dispersion_y = NAN,
				  const float corrected_dispersion_z = NAN)
			: particle_t(particle), _id(id), _tower_id(tower_id),
			  _partial_energy_sum(partial_energy_sum),
			  _electromagnetic_chi_square(electromagnetic_chi_square),
			  _tof_timing(tof_timing),
			  _incident_angle(incident_angle),
			  _corrected_dispersion_y(corrected_dispersion_y),
			  _corrected_dispersion_z(corrected_dispersion_z)
		{
		}
		cluster_t(const float energy, const float eta,
				  const float phi, const int id,
				  const std::vector<tower_id_t> tower_id =
				  std::vector<tower_id_t>(),
				  const std::vector<float> partial_energy_sum =
				  std::vector<float>(),
				  const float electromagnetic_chi_square = NAN,
				  const tof_timing_t tof_timing =
				  tof_timing_t(NAN, NAN, NAN, NAN),
				  const float incident_angle = NAN,
				  const float corrected_dispersion_y = NAN,
				  const float corrected_dispersion_z = NAN)
			: particle_t(energy, NAN, eta, phi), _id(id),
			  _tower_id(tower_id),
			  _partial_energy_sum(partial_energy_sum),
			  _electromagnetic_chi_square(electromagnetic_chi_square),
			  _tof_timing(tof_timing),
			  _incident_angle(incident_angle),
			  _corrected_dispersion_y(corrected_dispersion_y),
			  _corrected_dispersion_z(corrected_dispersion_z)
		{
		}
		cluster_t(const float momentum[4], const int id,
				  const std::vector<tower_id_t> tower_id,
				  const std::vector<float> partial_energy_sum,
				  const float electromagnetic_chi_square,
				  const float tof_timing[4],
				  const float incident_angle = NAN,
				  const float corrected_dispersion_y = NAN,
				  const float corrected_dispersion_z = NAN)
			: particle_t(momentum), _id(id), _tower_id(tower_id),
			  _partial_energy_sum(partial_energy_sum),
			  _electromagnetic_chi_square(electromagnetic_chi_square),
			  _tof_timing(tof_timing_t(tof_timing[0], tof_timing[1],
									   tof_timing[2], tof_timing[3])),
			  _incident_angle(incident_angle),
			  _corrected_dispersion_y(corrected_dispersion_y),
			  _corrected_dispersion_z(corrected_dispersion_z)
		{
		}
		inline int &id(void)
		{
			return _id;
		}
		inline int id(void) const
		{
			return _id;
		}
		inline std::vector<tower_id_t> &tower_id(void)
		{
			return _tower_id;
		}
		inline std::vector<tower_id_t> tower_id(void) const
		{
			return _tower_id;
		}
		inline std::vector<float> &partial_energy_sum(void)
		{
			return _partial_energy_sum;
		}
		inline std::vector<float> partial_energy_sum(void) const
		{
			return _partial_energy_sum;
		}
		inline tower_id_t central_tower_id(void) const
		{
			I(_tower_id.size() > 0);

			return _tower_id[0];
		}
		inline float &electromagnetic_chi_square(void)
		{
			return _electromagnetic_chi_square;
		}
		inline float electromagnetic_chi_square(void) const
		{
			return _electromagnetic_chi_square;
		}
		inline tof_timing_t &tof_timing(void)
		{
			return _tof_timing;
		}
		inline tof_timing_t tof_timing(void) const
		{
			return _tof_timing;
		}
		inline float &incident_angle(void)
		{
			return _incident_angle;
		}
		inline float incident_angle(void) const
		{
			return _incident_angle;
		}
		inline float &corrected_dispersion_y(void)
		{
			return _corrected_dispersion_y;
		}
		inline float corrected_dispersion_y(void) const
		{
			return _corrected_dispersion_y;
		}
		inline float &corrected_dispersion_z(void)
		{
			return _corrected_dispersion_z;
		}
		inline float corrected_dispersion_z(void) const
		{
			return _corrected_dispersion_z;
		}
		inline bool tower_id_match(const bool map[]) const
		{
			return map[central_tower_id()];
		}
		inline cluster_t energy_scale(const float scale[]) const
		{
			cluster_t retval = *this;

			retval.energy() *= scale[central_tower_id()];

			return retval;
		}
		friend class event_file_t;
	};

	class collisional_system_t;

	class nucleus_t {
	private:
		int _proton_number;
		int _nucleon_number;
	public:
		nucleus_t(void)
			: _proton_number(-1), _nucleon_number(-1)
		{
		}
		nucleus_t(const int proton_number, const int nucleon_number)
			: _proton_number(proton_number),
			  _nucleon_number(nucleon_number)
		{
			I(_proton_number >= 0);
			I(_nucleon_number >= 0);
			I(_nucleon_number >= _proton_number);
		}
		inline int proton_number(void) const
		{
			return _proton_number;
		}
		inline int nucleon_number(void) const
		{
			return _nucleon_number;
		}
		inline int neutron_number(void) const
		{
			return _nucleon_number - _proton_number;
		}
		inline bool is_hadron(void) const
		{
			I(_proton_number >= 0);
			I(_nucleon_number >= 0);
			I(_nucleon_number >= _proton_number);

			return _nucleon_number == 1;
		}
		std::string element_symbol(void) const;
		std::string isotope_symbol(void) const;
		inline std::string name(void) const
		{
			return isotope_symbol();
		}
		/**
		 * Returns the collisional (*this, target) system at rest
		 *
		 * @param[in] target the target nucleus
		 * @return the collisional (*this, target) system at rest
		 */
		collisional_system_t operator+(const nucleus_t target) const;
	};

	// Instances of RHIC and LHC projectiles

	// Both RHIC and LHC use/will use the proton
	extern nucleus_t proton;
	// RHIC uses the deuteron
	extern nucleus_t deuteron;
	// RHIC uses the Cu63 isotope, see: C. J. Gardner, "Booster, AGS,
	// and RHIC Parameters for the 2004-2005 RHIC Run" (2003).
	extern nucleus_t copper_63;
	// RHIC uses the Au197 isotope, which is also the only stable one
	extern nucleus_t gold_197;
	// Future LHC Pb208 isotope
	extern nucleus_t lead_208;
	// Future RHIC U238 isotope
	extern nucleus_t uranium_238;

	// The main RHIC operating energy, sqrt(s_NN) = 200 GeV
	extern double rhic_sqrt_s_nn;
	// The main LHC heavy ion operating energy, sqrt(s_NN) = 5520 GeV
	extern double lhc_sqrt_s_nn;
	// The main LHC hadronic operating energy, sqrt(s) = 1960 GeV
	extern double tevatron_sqrt_s;
	// The main LHC hadronic operating energy, sqrt(s) = 14000 GeV
	extern double lhc_sqrt_s;

	class nucleus_or_particle_t {
	private:
		nucleus_t _nucleus;
		particle_t _particle;
		bool _is_particle;
	};

	class collisional_system_t {
	private:
		// The "projectile/target" (versus Pythia's "beam/target")
		// naming follows the HIJING convention.
		nucleus_t _projectile;
		nucleus_t _target;
		float _cms_energy;
	public:
		collisional_system_t()
			: _projectile(), _target(), _cms_energy(0)
		{
			I(_cms_energy >= 0);
		}
		collisional_system_t(const nucleus_t projectile,
							 const nucleus_t target,
							 const float cms_energy)
			: _projectile(projectile), _target(target),
			  _cms_energy(cms_energy)
		{
			I(_cms_energy >= 0);
		}
		nucleus_t projectile(void) const
		{
			return _projectile;
		}
		nucleus_t target(void) const
		{
			return _target;
		}
		float cms_energy(void) const
		{
			return _cms_energy;
		}
		bool is_hadronic(void) const
		{
			return _projectile.is_hadron() && _target.is_hadron();
		}
		collisional_system_t at_cms_energy(const float cms_energy)
		{
			return collisional_system_t(_projectile, _target,
										cms_energy);
		}
		std::string name(void) const;
		friend class event_file_t;
	};

	inline collisional_system_t
	nucleus_t::operator+(const nucleus_t target) const
	{
		return collisional_system_t(*this, target, 0);
	}

	inline std::ostream &
	operator<<(std::ostream &os, const nucleus_t &nucleus)
	{
		os << nucleus.name();

		return os;
	}

	inline std::ostream &
	operator<<(std::ostream &os, const collisional_system_t &system)
	{
		os << system.name();

		return os;
	}

	template<typename physical_particle_list_t,
			 typename abstract_particle_t>
	inline physical_particle_list_t particle_filter(
		const physical_particle_list_t &particle_list,
		std::const_mem_fun_ref_t<bool, abstract_particle_t>
		member_function)
	{
		physical_particle_list_t retval;

		// Simple heuristics
		retval.reserve(particle_list.size() >> 2);

		for(typename physical_particle_list_t::const_iterator
				iterator = particle_list.begin();
			iterator != particle_list.end(); iterator++)
			if(member_function(*iterator))
				retval.push_back(*iterator);

		return retval;
	}

	template<typename physical_particle_list_t,
			 typename abstract_particle_t, typename property_t,
			 typename comparison_t>
	inline physical_particle_list_t particle_filter(
		const physical_particle_list_t &particle_list,
		std::binder2nd<comparison_t> comparison,
		std::const_mem_fun_ref_t<property_t, abstract_particle_t>
		member_function,
		const property_t reference)
	{
		physical_particle_list_t retval;

		// Simple heuristics
		retval.reserve(particle_list.size() >> 2);

		for(typename physical_particle_list_t::const_iterator
				iterator = particle_list.begin();
			iterator != particle_list.end(); iterator++)
			if(comparison(member_function(*iterator), reference))
				retval.push_back(*iterator);

		return retval;
	}

}

#endif // JETEVENT_PARTICLE_H_