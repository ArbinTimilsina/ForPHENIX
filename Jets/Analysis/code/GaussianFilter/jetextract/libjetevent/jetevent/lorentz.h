// -*- mode: c++; -*-

/////////////////////////////////////////////////////////////////////

#ifndef JETEVENT_LORENTZ_H_
#define JETEVENT_LORENTZ_H_

// LinkDef files generated by rootcint do not include config.h
#if defined(__CINT__) || defined(G__DICTIONARY)
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H
#endif // defined(__CINT__) || defined(G__DICTIONARY)

#include <jetbase/dbc.h>
#include <jetbase/num.h>
#include <iostream>

#if defined(HAVE_ROOT) && !defined(NVERIFY)
#include <TLorentzVector.h>
#endif // defined(HAVE_ROOT) && !defined(NVERIFY)

// JET RECONSTRUCTION NUMERICS - LORENTZ VECTOR KINEMATICS

namespace jet {

	class abstract_lorentz_vector_t {
	protected:
		static float
		max_pseudorapidity(const float perp, const float z);
		static float
		pseudorapidity_limit(const float perp,
							 const float pseudorapidity);
	public:
		virtual ~abstract_lorentz_vector_t(void)
		{
		}
		/**
		 * Returns the timelike or zeroth component of the Lorentz
		 * vector
		 *
		 * @return the timelike or zeroth component
		 * @see energy()
		 */
		virtual float time(void) const = 0;
		/**
		 * Returns x or the first Cartesian/spatial component of the
		 * Lorentz vector
		 *
		 * @return x or the first Cartesian/spatial component
		 */
		virtual float x(void) const = 0;
		/**
		 * Returns y or the second Cartesian/spatial component of the
		 * Lorentz vector
		 *
		 * @return y or the second Cartesian/spatial component
		 */
		virtual float y(void) const = 0;
		/**
		 * Returns z or the third Cartesian/spatial component of the
		 * Lorentz vector
		 *
		 * @return z or the third Cartesian/spatial component
		 */
		virtual float z(void) const = 0;
		/**
		 * Returns the timelike or zeroth component of the Lorentz
		 * vector (alias for time())
		 *
		 * @return the timelike or zeroth component
		 * @see time()
		 */
		inline float energy(void) const
		{
			return time();
		}
		/**
		 * Returns the square of the transverse component (with
		 * respect to z)
		 *
		 * @return the square of the transverse component
		 * @see perp()
		 */
		virtual float perp_square(void) const;
		/**
		 * Returns the transverse component (with respect to z)
		 *
		 * @return the transverse component
		 * @see perp_square()
		 */
		virtual float perp(void) const;
		// Note: cartesian_magnitude() is not called euclidean_norm()
		// to avoid confusing with a 4D Euclidean norm
		/**
		 * Returns the square of the magnitude/Euclidean norm of the
		 * Cartesian/spatial components
		 *
		 * @return the square of the magnitude/Euclidean norm of the
		 * Cartesian/spatial components
		 * @see cartesian_magnitude()
		 */
		virtual float cartesian_magnitude_square(void) const;
		/**
		 * Returns the magnitude/Euclidean norm of the
		 * Cartesian/spatial components
		 *
		 * @return the magnitude/Euclidean norm of the
		 * Cartesian/spatial components
		 * @see cartesian_magnitude_square()
		 */
		virtual float cartesian_magnitude(void) const;
		// Note: The method norm() is intentially avoided to prevent
		// ambiguities
		// FIXME: Rename this to minkowski_norm() and make mass() an
		// alias
		/**
		 * Returns the square of the mass/Minkowski norm
		 *
		 * @return the square of the mass/Minkowski norm
		 * @see mass()
		 */
		virtual float mass_square(void) const;
		/**
		 * Returns the mass/Minkowski norm
		 *
		 * @return the mass/Minkowski norm
		 * @see mass_square()
		 */
		virtual float mass(void) const;
		virtual float mass_perp_square(void) const;
		virtual float mass_perp(void) const;
		virtual float pseudorapidity(void) const;
		virtual float rapidity(void) const;
		virtual float sin_polar_angle(void) const;
		virtual float time_perp(void) const;
		inline float energy_perp(void) const
		{
			return time_perp();
		}
		virtual float azimuth(void) const;
		virtual float
		cartesian_dot(const abstract_lorentz_vector_t &v) const;
		template<typename vector_t>
		inline float longitudinal_magnitude(const vector_t &v) const
		{
			return cartesian_dot(v) / v.cartesian_magnitude();
		}
		template<typename vector_t>
		inline float longitudinal_fraction(const vector_t &v) const
		{
			return cartesian_dot(v) / v.cartesian_magnitude_square();
		}
		template<typename vector_t>
		inline float cos_angle(const vector_t &v) const
		{
			return cartesian_dot(v) /
				(cartesian_magnitude() * v.cartesian_magnitude());
		}
		float pseudorapidity_edge_distance(
			const float pseudorapidity_min,
			const float pseudorapidity_max) const;
		float azimuth_edge_distance(
			const float azimuth_min,
			const float azimuth_max) const;
		float azimuth_edge_distance(
			const float azimuth_min_0,
			const float azimuth_max_0,
			const float azimuth_min_1,
			const float azimuth_max_1) const;
	};

    class event_file_t;

#ifdef HAVE_ROOT
	// Sadly, CINT does not support any in-language alignment
	// specification, i.e. neither the GNU style __attribute__, nor
	// the Microsoft/Intel __declspec.
#ifdef __INTEL_COMPILER
#pragma pack(push)
#endif // __INTEL_COMPILER
#pragma pack(16)
#endif // HAVE_ROOT

#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 810)
#pragma warning(disable: 1572)
#endif // __INTEL_COMPILER
    class lorentz_vector_t : public abstract_lorentz_vector_t {
    protected:
		float _x[4];
    public:
		lorentz_vector_t(void)
		{
		}
		lorentz_vector_t(const float time, const float x,
						 const float y, const float z);
		lorentz_vector_t(const float time, const float x[3]);

		inline void set_energy(const float e)
		{
			I(_x != NULL);
#ifndef HAVE_PHENIX
			I(F(e));
#endif // HAVE_PHENIX
			I(e >= 0);

			_x[0] = e;
		}

		inline float operator[](const int n) const
		{
			return _x[n];
		}

		inline float &operator[](const int n)
		{
			return _x[n];
		}

		inline float time(void) const
		{
			return _x[0];
		}

		inline float &time(void)
		{
			return _x[0];
		}

		inline float &energy(void)
		{
			return _x[0];
		}

		inline float x(void) const
		{
			return _x[1];
		}

		inline float &x(void)
		{
			return _x[1];
		}

		inline float y(void) const
		{
			return _x[2];
		}

		inline float &y(void)
		{
			return _x[2];
		}

		inline float z(void) const
		{
			return _x[3];
		}

		inline float &z(void)
		{
			return _x[3];
		}

		inline lorentz_vector_t operator*(const float s) const
		{
			return lorentz_vector_t(_x[0] * s, _x[1] * s,
									_x[2] * s, _x[3] * s);
		}

		inline float perp(void) const
		{
			return sqrtf(_x[1] * _x[1] + _x[2] * _x[2]);
		}

		inline float azimuth(void) const
		{
			return atan2f(_x[2], _x[1]);
		}

		inline float height(void) const
		{
			return z();
		}

		inline lorentz_vector_t
		operator+(const lorentz_vector_t &v) const
		{
			return lorentz_vector_t(
				_x[0] + v._x[0], _x[1] + v._x[1],
				_x[2] + v._x[2], _x[3] + v._x[3]);
		}

		inline lorentz_vector_t
		operator+=(const lorentz_vector_t &v)
		{
			_x[0] += v._x[0];
			_x[1] += v._x[1];
			_x[2] += v._x[2];
			_x[3] += v._x[3];

			return *this;
		}

		inline lorentz_vector_t
		operator-(const lorentz_vector_t &v) const
		{
			return lorentz_vector_t(
				_x[0] - v._x[0], _x[1] - v._x[1],
				_x[2] - v._x[2], _x[3] - v._x[3]);
		}

		inline lorentz_vector_t
		operator-=(const lorentz_vector_t &v)
		{
			_x[0] -= v._x[0];
			_x[1] -= v._x[1];
			_x[2] -= v._x[2];
			_x[3] -= v._x[3];

			return *this;
		}

		lorentz_vector_t transverse(const lorentz_vector_t &v) const;
		float transverse_magnitude(const lorentz_vector_t &v) const;

		inline bool operator==(const float c) const
		{
			return (_x[0] == c && _x[1] == c &&
					_x[2] == c && _x[3] == c);
		}

		inline bool operator!=(const float c) const
		{
			return (_x[0] != c || _x[1] != c ||
					_x[2] != c || _x[3] != c);
		}

		friend lorentz_vector_t
		operator*(const float s, const lorentz_vector_t v);
		friend class event_file_t;
#ifdef __CINT__
	};
#else // __CINT__
	}
#ifndef HAVE_ROOT
	__attribute__ ((aligned(16)))
#endif // HAVE_ROOT
	;
#endif // __CINT__
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif // __INTEL_COMPILER

	inline lorentz_vector_t
	operator*(const float s, const lorentz_vector_t v)
	{
		return v * s;
	}

#ifdef HAVE_ROOT
#ifdef __INTEL_COMPILER
#pragma pack(pop)
#endif // __INTEL_COMPILER
#endif // HAVE_ROOT

}

#endif // JETEVENT_LORENTZ_H_