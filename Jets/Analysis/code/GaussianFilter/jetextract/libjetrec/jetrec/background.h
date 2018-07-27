// -*- mode: c++; -*-

#ifndef XJETREC_BACKGROUND_H_
#define XJETREC_BACKGROUND_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <jetrec/series.h>

/////////////////////////////////////////////////////////////////////

// Background Model

namespace jet {

#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 869)
#endif // __INTEL_COMPILER

	class collision_geometry_t {
	private:
		float _vertex_height;
		float _centrality;
		float _reaction_plane;
	public:
		inline collision_geometry_t(void)
			: _vertex_height(NAN), _centrality(NAN),
			  _reaction_plane(NAN)
		{
		}
		inline collision_geometry_t(
			const float vertex_height, const float centrality,
			const float reaction_plane)
			: _vertex_height(vertex_height),
			  _centrality(centrality),
			  _reaction_plane(reaction_plane)
		{
		}
		inline float vertex_height(void) const
		{
			return _vertex_height;
		}
		inline float &vertex_height(void)
		{
			return _vertex_height;
		}
		inline float centrality(void) const
		{
			return _centrality;
		}
		inline float &centrality(void)
		{
			return _centrality;
		}
		inline float reaction_plane(void) const
		{
			return _reaction_plane;
		}
		inline float &reaction_plane(void)
		{
			return _reaction_plane;
		}
		inline bool in_range(void) const
		{
			const bool ret =
				_vertex_height >= -25.0F && _vertex_height < 25.0F &&
				_centrality >= 0.0F && _centrality <= 100.0F;

#if 0
			if(!ret) {
				std::cerr << __FILE__ << ':' << __LINE__
						  << ": _vertex_height = " << _vertex_height
						  << ", _centrality = " << _centrality
						  << std::endl;
			}
#endif

			return ret;
		}
	};

	class background_model_t {
	public:
		static const unsigned int nbin_centrality = 101;
		static const unsigned int nbin_vertex = 10;
		static const unsigned int nbin_reaction_plane = 32;
	protected:
		bool _factorized;
		inline background_model_t(const bool factorized)
			: _factorized(factorized)
		{
		}
		inline virtual ~background_model_t(void)
		{
		}
	public:
		inline bool factorized(void) const
		{
			return _factorized;
		}
		inline virtual double operator()(
			const collision_geometry_t &geometry,
			const double pseudorapidity, const double azimuth) const
		{
			return 0;
		}
		inline virtual void evaluate_order_2(
			double &f, double gradient[], double hessian[],
			const collision_geometry_t &geometry,
			const double pseudorapidity, const double azimuth) const
		{
		}
		inline virtual const float *
		vertex_centrality_dependence(void) const
		{
			return NULL;
		}
		inline virtual const float vertex_centrality_dependence(
			const collision_geometry_t &geometry) const
		{
			return 0;
		}
		inline virtual const chebyshev_series_nd_sparse_t *
		position_dependence(void) const
		{
			return NULL;
		}
		inline virtual chebyshev_series_nd_sparse_t
		position_dependence(
			const collision_geometry_t &geometry) const
		{
			return chebyshev_series_nd_sparse_t();
		}
	};

#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 1125)
#endif // __INTEL_COMPILER

	class factorized_background_model_t : public background_model_t {
	public:
		static const unsigned int nbin_centrality = 101;
		static const unsigned int nbin_vertex = 10;
	private:
		float _vertex_centrality_dependence[
			nbin_vertex * nbin_centrality];
		chebyshev_series_nd_sparse_t
		_position_dependence[nbin_vertex];
	public:
		static unsigned int vertex_discretize(
			const collision_geometry_t &geometry);
		static unsigned int vertex_centrality_discretize(
			const collision_geometry_t &geometry);
		inline factorized_background_model_t(void)
			: background_model_t(true)
		{
			for(unsigned int i = 0;
				i < nbin_vertex * nbin_centrality; i++) {
				_vertex_centrality_dependence[i] = NAN;
			}
		}
		inline factorized_background_model_t(
			const float vertex_centrality_dependence[],
			const chebyshev_series_nd_sparse_t position_dependence[])
			: background_model_t(true)
		{
			memcpy(_vertex_centrality_dependence,
				   vertex_centrality_dependence,
				   nbin_vertex * nbin_centrality * sizeof(float));
			for(unsigned int i = 0; i < nbin_vertex; i++) {
				_position_dependence[i] = position_dependence[i];
			}
		}
		inline const float *vertex_centrality_dependence(void) const
		{
			return _vertex_centrality_dependence;
		}
		inline float *vertex_centrality_dependence(void)
		{
			return _vertex_centrality_dependence;
		}
		inline const chebyshev_series_nd_sparse_t *
		position_dependence(void) const
		{
			return _position_dependence;
		}
		inline chebyshev_series_nd_sparse_t *
		position_dependence(void)
		{
			return _position_dependence;
		}
		inline chebyshev_series_nd_sparse_t
		position_dependence(
			const unsigned int vertex,
			const unsigned int centrality,
			const unsigned int reaction_plane) const
		{
			return _position_dependence[vertex];
		}
		inline const float vertex_centrality_dependence(
			const collision_geometry_t &geometry) const
		{
			return _vertex_centrality_dependence[
				vertex_centrality_discretize(geometry)];;
		}
		inline chebyshev_series_nd_sparse_t
		position_dependence(
			const collision_geometry_t &geometry) const
		{
			return _position_dependence[
				vertex_discretize(geometry)];
		}
		inline bool empty(void) const
		{
			return _position_dependence[0].empty();
		}
		double operator()(
			const collision_geometry_t &geometry,
			const double pseudorapidity, const double azimuth) const;
		void evaluate_order_2(
			double &f, double gradient[], double hessian[],
			const collision_geometry_t &geometry,
			const double pseudorapidity, const double azimuth) const;
	};

	class general_background_model_t : public background_model_t {
	private:
		chebyshev_series_nd_sparse_t _position_dependence[
			nbin_vertex * nbin_centrality * nbin_reaction_plane];
	public:
		static void discretize_reaction_plane(
			unsigned int &index_low, unsigned int &index_high,
			double &weight_low, double &weight_high,
			const double reaction_plane);
		inline general_background_model_t(void)
			: background_model_t(false)
		{
		}
		inline general_background_model_t(
			const chebyshev_series_nd_sparse_t position_dependence[])
			: background_model_t(false)
		{
			for(unsigned int i = 0;
				i < nbin_centrality * nbin_reaction_plane; i++) {
				_position_dependence[i] = position_dependence[i];
			}
		}
		inline virtual const chebyshev_series_nd_sparse_t *
		position_dependence(void) const
		{
			return _position_dependence;
		}
		inline chebyshev_series_nd_sparse_t
		position_dependence(
			const unsigned int vertex,
			const unsigned int centrality,
			const unsigned int reaction_plane) const
		{
			return _position_dependence[
				(vertex * nbin_centrality + centrality) *
				nbin_reaction_plane + reaction_plane];
		}
		inline chebyshev_series_nd_sparse_t &
		position_dependence(
			const unsigned int vertex,
			const unsigned int centrality,
			const unsigned int reaction_plane)
		{
			return _position_dependence[
				(vertex * nbin_centrality + centrality) *
				nbin_reaction_plane + reaction_plane];
		}
		inline chebyshev_series_nd_sparse_t
		position_dependence(
			const collision_geometry_t &geometry) const
		{
			// FIXME: not implemented
			return chebyshev_series_nd_sparse_t();
		}
		double operator()(
			const double vertex, const double centrality,
			const double reaction_plane, const double pseudorapidity,
			const double azimuth) const;
		void evaluate_order_2(
			double &f, double gradient[], double hessian[],
			const double vertex, const double centrality,
			const double reaction_plane, const double pseudorapidity,
			const double azimuth) const;
	};

#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif // __INTEL_COMPILER

#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif // __INTEL_COMPILER

}

#endif // XJETREC_BACKGROUND_H_
