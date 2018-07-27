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

#ifdef HAVE_KTJET
	KtJet::KtLorentzVector reconstruction_ktjet_t::
	kt_lorentz_vector(const track_t &track)
	{
		return KtJet::KtLorentzVector(
			track.momentum().x(),
			track.momentum().y(),
			track.momentum().z(),
			track.energy());
	}

	snowmass_vector_t reconstruction_ktjet_t::
	snowmass_vector(
		const KtJet::KtLorentzVector &kt_lorentz_vector)
	{
		jet::snowmass_vector_t momentum;

		momentum.set_time((float)kt_lorentz_vector.e());
		momentum.set_x_y_z((float)kt_lorentz_vector.px(),
						   (float)kt_lorentz_vector.py(),
						   (float)kt_lorentz_vector.pz());

		return momentum;
	}

	std::vector<jet_t> reconstruction_ktjet_t::
	reconstruct(
		const std::vector<track_t>::const_iterator &track_begin,
		const std::vector<track_t>::const_iterator &track_end,
		const bool refine)
	{
		std::vector<KtJet::KtLorentzVector> pseudojet;

		for(std::vector<track_t>::const_iterator iterator =
				track_begin;
			iterator != track_end; iterator++)
#if 0
			if(iterator->is_final_state())
#endif
				pseudojet.push_back(kt_lorentz_vector(*iterator));

		KtJet::KtEvent event(pseudojet, _collision_type,
							 _distance_scheme, _recombination_scheme,
							 _radius);
		std::vector<KtJet::KtLorentzVector> jet = event.getJetsPt();
		std::vector<jet_t> ret;

		for(std::vector<KtJet::KtLorentzVector>::const_iterator
				iterator = jet.begin();
			iterator != jet.end(); iterator++) {
			if(iterator->perp() > FLT_EPSILON)
				ret.push_back(snowmass_vector(*iterator));
		}

		return ret;
	}
#endif // HAVE_KTJET

#ifdef HAVE_FASTJET
	fastjet::PseudoJet reconstruction_fastjet_t::
	pseudo_jet(const track_t &track)
	{
		if(track.momentum().perp() < FLT_MIN) {
			snowmass_vector_t modified_momentum = track.momentum();

			modified_momentum.perp() = FLT_MIN;
			return fastjet::PseudoJet(
				modified_momentum.x(), modified_momentum.y(),
				modified_momentum.z(), modified_momentum.time());
		}
		else
			return fastjet::PseudoJet(
				track.momentum().x(), track.momentum().y(),
				track.momentum().z(), track.energy());
	}

	snowmass_vector_t reconstruction_fastjet_t::
	snowmass_vector(const fastjet::PseudoJet &pseudo_jet)
	{
		jet::snowmass_vector_t momentum;

		momentum.set_time((float)pseudo_jet.e());
		momentum.set_x_y_z((float)pseudo_jet.px(),
						   (float)pseudo_jet.py(),
						   (float)pseudo_jet.pz());

		return momentum;
	}

	std::vector<jet_t> reconstruction_fastjet_t::
	reconstruct(
		const std::vector<track_t>::const_iterator &track_begin,
		const std::vector<track_t>::const_iterator &track_end,
		const bool refine)
	{
		if(track_begin == track_end) {
			return std::vector<jet_t>();
		}

		std::vector<fastjet::PseudoJet> pseudojet;

		for(std::vector<track_t>::const_iterator iterator =
				track_begin;
			iterator != track_end; iterator++) {
			fastjet::PseudoJet p = pseudo_jet(*iterator);

			p.set_user_index(iterator - track_begin);
			pseudojet.push_back(p);
		}

		fastjet::ClusterSequence
			cluster_sequence(pseudojet, _jet_definition);
		std::vector<fastjet::PseudoJet> jet =
			cluster_sequence.inclusive_jets(0);

#if 0
		for(std::vector<fastjet::PseudoJet>::const_iterator
				iterator_jet = jet.begin();
			iterator_jet != jet.end(); iterator_jet++) {
			const jet_t p = snowmass_vector(*iterator_jet);
			std::vector<fastjet::PseudoJet> constituent =
				cluster_sequence.constituents(*iterator_jet);

			for(std::vector<fastjet::PseudoJet>::const_iterator
					iterator_constituent = constituent.begin();
				iterator_constituent != constituent.end();
				iterator_constituent++) {
				std::vector<track_t>::const_iterator t =
					track_begin +
					iterator_constituent->user_index();

				if(t->momentum().radial_distance(p.momentum()) <
				   _resummation_radius) {
					fprintf(stderr, "<FJT> %.8e %.8e %ld\n",
							iterator_constituent->eta(),
							jet::angular_range_reduce(iterator_constituent->phi()),
							iterator_jet - jet.begin());
				}
			}
		}
#endif

		// Recover the heavy ion energy via resummantion (e.g.
		// SISCone)

		std::vector<jet_t> ret;

		for(std::vector<fastjet::PseudoJet>::const_iterator
				iterator = jet.begin();
			iterator != jet.end(); iterator++) {
			if(iterator->perp() > FLT_EPSILON) {
				const std::vector<fastjet::PseudoJet> constituent =
					cluster_sequence.constituents(*iterator);
				jet_t p = snowmass_vector(*iterator);

				p.constituent().reserve(constituent.size());

				float perp = 0;

				for(std::vector<fastjet::PseudoJet>::const_iterator
						constituent_iterator = constituent.begin();
					constituent_iterator != constituent.end();
					constituent_iterator++) {
					std::vector<track_t>::const_iterator t =
						track_begin +
						constituent_iterator->user_index();

					if(t->momentum().radial_distance(p.momentum()) <
					   _resummation_radius) {
						p.constituent().push_back(t);
						perp += t->momentum().perp();
					}
				}
				p.momentum().perp() = perp;
				ret.push_back(p);
			}
		}

		return ret;
	}

	double reconstruction_fastkt_area_t::
	perp_per_area_scale(const double x) const
	{
		const double f =
89.89811808379382 + x*(-3.7063608200282303 + 
      x*(3.063302726514298 + 
         x*(1.1614171035669132 + 
            x*(2.2267769146326164 + 
               x*(-0.43505463489965657 + 
                  x*(-2.409901965705066 + 
                     x*
                     (0.13929969004127662 + 
                     x*
                     (0.6619024158364079 + 
                     x*
                     (-0.027547577630292364 + 
                     x*
                     (-0.09861065034810801 + 
                     x*
                     (0.00344386607440872 + 
                     x*
                     (0.009400049449826028 + 
                     x*
                     (-0.0002884707537046019 + 
                     x*
                     (-0.0006136347515912258 + 
                     x*
                     (0.0000168833809566842 + 
                     x*
                     (0.000028376852894482593 + 
                     x*
                     (-7.085018622335271e-7 + 
                     x*
                     (-9.447690290827544e-7 + 
                     x*
                     (2.159786539854795e-8 + 
                     x*
                     (2.272208565603866e-8 + 
                     x*
                     (-4.790852487046989e-10 + 
                     x*
                     (-3.9110167139213594e-10 + 
                     x*
                     (7.653144334137503e-12 + 
                     x*
                     (4.696397344557738e-12 + 
                     x*
                     (-8.57547103793019e-14 + 
                     x*
                     (-3.735478067509467e-14 + 
                     x*
                     (6.395284659176343e-16 + 
                     x*
                     (1.768112504257168e-16 + 
                     x*
                     (-2.8503000363095586e-18 + 
                     (-3.769316720408583e-19 + 
                     5.7433132407915245e-21*x)*x)))))))))
                     ))))))))))))))))))));

		// Normalize by integral between [-3, 3]
		return f * (1.0 / 81.85115255685353);
	}

	double reconstruction_fastkt_area_t::
	sigma_perp_per_area_scale(const double x) const
	{
		const double f =
26.57796901859253 + x*(2.2189390214363747 + 
      x*(-1.5323018054755466 + 
         x*(-0.04062097510428009 + 
            x*(-1.0884609912388754 + 
               x*(-0.1723197240494518 + 
                  x*(0.7347816026228225 + 
                     x*
                     (0.05588915519493745 + 
                     x*
                     (-0.220611664977097 + 
                     x*
                     (-0.00998052966553496 + 
                     x*
                     (0.037208063722026286 + 
                     x*
                     (0.001171024257429219 + 
                     x*
                     (-0.003928873186745211 + 
                     x*
                     (-0.00009503650898606608 + 
                     x*
                     (0.0002772690086748692 + 
                     x*
                     (5.48181000955765e-6 + 
                     x*
                     (-0.000013599861141428215 + 
                     x*
                     (-2.2860480209716468e-7 + 
                     x*
                     (4.7360436585157485e-7 + 
                     x*
                     (6.952036039325521e-9 + 
                     x*
                     (-1.1793679219786107e-8 + 
                     x*
                     (-1.540991543530081e-10 + 
                     x*
                     (2.0862025158788238e-10 + 
                     x*
                     (2.4613586101922747e-12 + 
                     x*
                     (-2.560169904040913e-12 + 
                     x*
                     (-2.757835491918989e-14 + 
                     x*
                     (2.072239210054306e-14 + 
                     x*
                     (2.0561640572893647e-16 + 
                     x*
                     (-9.948709976174868e-17 + 
                     x*
                     (-9.158686871880626e-19 + 
                     (2.14568787134882e-19 + 
                     1.843675927767788e-21*x)*x))))))))))
                     )))))))))))))))))));

		// Normalize by integral between [-3, 3]
		return f * (1.0 / 21.098678631663283);
	}

	std::vector<jet_t> reconstruction_fastkt_area_t::
	reconstruct(
		const std::vector<track_t>::const_iterator &track_begin,
		const std::vector<track_t>::const_iterator &track_end,
		const bool refine)
	{
		if(track_begin == track_end) {
			return std::vector<jet_t>();
		}

		std::vector<fastjet::PseudoJet> pseudojet;

		for(std::vector<track_t>::const_iterator iterator =
				track_begin;
			iterator != track_end; iterator++) {
			fastjet::PseudoJet p = pseudo_jet(*iterator);

			p.set_user_index(iterator - track_begin);
			pseudojet.push_back(p);
		}

		const fastjet::ClusterSequenceArea
			cluster_sequence(pseudojet, _jet_definition,
							 _area_definition);
		const std::vector<fastjet::PseudoJet> jet =
			cluster_sequence.inclusive_jets(0);
		std::vector<fastjet::PseudoJet> fake_jet = jet;
		long previous_fake_jet_size = -1;
		double mean_perp_per_area = NAN;
		double standard_deviation_perp_per_area = NAN;

#if 0
		for(std::vector<fastjet::PseudoJet>::const_iterator
				iterator_jet = jet.begin();
			iterator_jet != jet.end(); iterator_jet++) {
			std::vector<fastjet::PseudoJet> constituent =
				cluster_sequence.constituents(*iterator_jet);

			for(std::vector<fastjet::PseudoJet>::const_iterator
					iterator_constituent = constituent.begin();
				iterator_constituent != constituent.end();
				iterator_constituent++) {
				fprintf(stderr, "<FJT> %.8e %.8e %ld\n",
						iterator_constituent->eta(),
						jet::angular_range_reduce(iterator_constituent->phi()),
						iterator_jet - jet.begin());
			}
		}
#endif

		while((long)fake_jet.size() != previous_fake_jet_size) {
			double sum_perp_per_area = 0;
			double sum_perp_per_area_square = 0;
			unsigned long count = 0;

			for(std::vector<fastjet::PseudoJet>::const_iterator
					iterator = fake_jet.begin();
				iterator != fake_jet.end(); iterator++) {
				if(fabs(iterator->pseudorapidity()) <= 3.0) {
					const double perp_per_area = iterator->perp() /
						cluster_sequence.area(*iterator);

					if(std::isfinite(perp_per_area)) {
						sum_perp_per_area += perp_per_area;
						sum_perp_per_area_square +=
							perp_per_area * perp_per_area;
						count++;
					}
				}
			}
			mean_perp_per_area = sum_perp_per_area / count;
			standard_deviation_perp_per_area = sqrt(std::max(
				0.0,
				sum_perp_per_area_square / count -
				mean_perp_per_area * mean_perp_per_area));

			std::vector<fastjet::PseudoJet> new_fake_jet;

			for(std::vector<fastjet::PseudoJet>::const_iterator
					iterator = fake_jet.begin();
				iterator != fake_jet.end(); iterator++) {
				const double area = cluster_sequence.area(*iterator);
				const double sigma_perp_per_area =
					(iterator->perp() -
					 area * mean_perp_per_area *
					 perp_per_area_scale(
						iterator->pseudorapidity())) /
					(area * standard_deviation_perp_per_area *
					 sigma_perp_per_area_scale(
						iterator->pseudorapidity()));

				if(sigma_perp_per_area <= 3) {
					new_fake_jet.push_back(*iterator);
				}
			}
#if 0
			fprintf(stderr, "%s:%d: fake_jet.size() = %lu%s\n",
					__FILE__, __LINE__, fake_jet.size(),
					new_fake_jet.size() == fake_jet.size() ?
					", terminating" : "");
#endif
			previous_fake_jet_size = fake_jet.size();
			fake_jet = new_fake_jet;
			if(!_background_subtract) {
				break;
			}
			// FIXME: I need to add a switch for fake rate study, in
			// which case there should be no iteration over fake jets!
		}

		// Recover the heavy ion energy via resummantion (e.g.
		// SISCone)

		std::vector<jet_t> ret;

		for(std::vector<fastjet::PseudoJet>::const_iterator
				iterator = jet.begin();
			iterator != jet.end(); iterator++) {
			if(iterator->perp() > FLT_EPSILON) {
				const std::vector<fastjet::PseudoJet> constituent =
					cluster_sequence.constituents(*iterator);
				jet_t p = snowmass_vector(*iterator);

				p.constituent().reserve(constituent.size());
				for(std::vector<fastjet::PseudoJet>::const_iterator
						constituent_iterator = constituent.begin();
					constituent_iterator != constituent.end();
					constituent_iterator++) {
					std::vector<track_t>::const_iterator t =
						track_begin +
						constituent_iterator->user_index();

					if(t->momentum().radial_distance(p.momentum()) <
					   _resummation_radius) {
						p.constituent().push_back(t);
					}
				}
				p.area() = cluster_sequence.area(*iterator);
#if 0
				fprintf(stderr, "%s:%d: %.8e %.8e\n",
						__FILE__, __LINE__,
						p.momentum().pseudorapidity(),
						p.momentum().perp() / p.area());
#endif
				if(_background_subtract) {
					p.momentum().perp() -=
						p.area() * mean_perp_per_area *
						perp_per_area_scale(
							p.momentum().pseudorapidity());
				}
				p.sigma_perp_per_area() = p.momentum().perp() /
					(p.area() * standard_deviation_perp_per_area *
					 sigma_perp_per_area_scale(
						p.momentum().pseudorapidity()));
				ret.push_back(p);
			}
		}

		return ret;
	}
#endif // HAVE_FASTJET

}
