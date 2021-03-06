// -*- mode: c++; -*-

#ifndef JETBASE_POLAR_H_
#define JETBASE_POLAR_H_

/////////////////////////////////////////////////////////////////////

// LinkDef files generated by rootcint do not include config.h
#if defined(__CINT__) || defined(G__DICTIONARY)
#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H
#endif // defined(__CINT__) || defined(G__DICTIONARY)

#include <jetbase/dbc.h>
#include <jetbase/specfunc.h>

namespace jet {

#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 810)
#pragma warning(disable: 1572)
#endif // __INTEL_COMPILER
	class cartesian_point_t {
	private:
		float _x[2];
	public:
		inline cartesian_point_t(void)
		{
		}
		inline cartesian_point_t(const cartesian_point_t &point)
		{
			_x[0] = point._x[0];
			_x[1] = point._x[1];
		}
		inline cartesian_point_t(const float x, const float y)
		{
			_x[0] = x;
			_x[1] = y;
		}
		inline cartesian_point_t(const float x[2])
		{
			_x[0] = x[0];
			_x[1] = x[1];
		}
		inline const float *x(void) const
		{
			return _x;
		}
		inline float *x(void)
		{
			return _x;
		}
		inline float operator[](const int n) const
		{
			return _x[n];
		}
		inline float &operator[](const int n)
		{
			return _x[n];
		}
		inline cartesian_point_t
		operator+(const cartesian_point_t &point) const
		{
			return cartesian_point_t(_x[0] + point._x[0],
									 _x[1] + point._x[1]);
		}
		inline cartesian_point_t
		operator-(const cartesian_point_t &point) const
		{
			return cartesian_point_t(_x[0] - point._x[0],
									 _x[1] - point._x[1]);
		}
		inline cartesian_point_t
		operator+=(const cartesian_point_t &point)
		{
			_x[0] += point._x[0];
			_x[1] += point._x[1];

			return *this;
		}
		inline cartesian_point_t
		operator-=(const cartesian_point_t &point)
		{
			_x[0] -= point._x[0];
			_x[1] -= point._x[1];

			return *this;
		}
		inline cartesian_point_t operator*(const float scale) const
		{
			return cartesian_point_t(_x[0] * scale, _x[1] * scale);
		}
		inline cartesian_point_t operator/(const float scale) const
		{
			return cartesian_point_t(_x[0] / scale, _x[1] / scale);
		}
		inline float dot(const cartesian_point_t &point) const
		{
			return _x[0] * point._x[0] + _x[1] * point._x[1];
		}
		inline float cross(const cartesian_point_t &point) const
		{
			return _x[0] * point._x[1] - _x[1] * point._x[0];
		}
		inline float norm_square(void) const
		{
			return _x[0] * _x[0] + _x[1] * _x[1];
		}
		inline float norm(void) const
		{
			return sqrtf(norm_square());
		}
		inline cartesian_point_t unit_point(void) const
		{
			return *this / norm();
		}
		inline cartesian_point_t rotate_cw(void) const
		{
			return cartesian_point_t(_x[1], -_x[0]);
		}
		inline cartesian_point_t rotate_ccw(void) const
		{
			return cartesian_point_t(-_x[1], _x[0]);
		}
		inline bool operator==(const cartesian_point_t &point) const
		{
			return _x[0] == point._x[0] && _x[1] == point._x[1];
		}
		inline bool operator!=(const cartesian_point_t &point) const
		{
			return _x[0] != point._x[0] || _x[1] != point._x[1];
		}
		friend cartesian_point_t
		operator*(const float scale, const cartesian_point_t &point);
	};
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif // __INTEL_COMPILER

	inline cartesian_point_t operator*(const float scale,
									   const cartesian_point_t &point)
	{
		return cartesian_point_t(scale * point._x[0],
								 scale * point._x[1]);
	}

	class cylindrical_point_t {
    protected:
		float _azimuth;
		float _height;
    public:
		cylindrical_point_t(void)
		{
		}
		cylindrical_point_t(const float azimuth, const float height)
			: _azimuth(azimuth), _height(height)
		{
		}
		cylindrical_point_t(const float x[2])
			: _azimuth(x[0]), _height(x[1])
		{
		}
		inline float &azimuth(void)
		{
			return _azimuth;
		}
		inline float azimuth(void) const
		{
			return _azimuth;
		}
		inline float &height(void)
		{
			return _height;
		}
		inline float height(void) const
		{
			return _height;
		}
		inline float radius(void) const
		{
			return sqrtf(_azimuth * _azimuth + _height * _height);
		}
    };

	class angular_slope_t {
    protected:
		float _azimuth;
		float _inclination;
    public:
		angular_slope_t(void)
		{
		}
		angular_slope_t(const float azimuth, const float inclination)
			: _azimuth(azimuth), _inclination(inclination)
		{
		}
		angular_slope_t(const float x[2])
			: _azimuth(x[0]), _inclination(x[1])
		{
		}
		inline float &azimuth(void)
		{
			return _azimuth;
		}
		inline float azimuth(void) const
		{
			return _azimuth;
		}
		inline float &inclination(void)
		{
			return _inclination;
		}
		inline float inclination(void) const
		{
			return _inclination;
		}
    };

    class spherical_point_t {
    protected:
		float _polar_angle;
		float _azimuth;
    public:
		spherical_point_t(void)
		{
		}
		spherical_point_t(const float polar_angle,
						  const float azimuth)
			: _polar_angle(polar_angle), _azimuth(azimuth)
		{
			IG(_polar_angle >= 0 && _polar_angle <= M_PI, F(_polar_angle));
			IG(_azimuth >= -M_PI && _azimuth <= M_PI, F(_azimuth));
		}
		spherical_point_t(const float x[2])
			: _polar_angle(x[0]), _azimuth(x[1])
		{
			IG(_polar_angle >= 0 && _polar_angle <= M_PI, F(_polar_angle));
			IG(_azimuth >= -M_PI && _azimuth <= M_PI, F(_azimuth));
		}
		inline float &polar_angle(void)
		{
			return _polar_angle;
		}
		inline float polar_angle(void) const
		{
			I(_polar_angle >= 0 && _polar_angle <= M_PI);

			return _polar_angle;
		}
		inline float &azimuth(void)
		{
			return _azimuth;
		}
		inline float azimuth(void) const
		{
			I(_azimuth >= -M_PI && _azimuth <= M_PI);

			return _azimuth;
		}
		inline float central_angle(const spherical_point_t &p) const
		{
			I(_polar_angle >= 0 && _polar_angle <= M_PI);
			I(_azimuth >= -M_PI && _azimuth <= M_PI);
			I(p._polar_angle >= 0 && p._polar_angle <= M_PI);
			I(p._azimuth >= -M_PI && p._azimuth <= M_PI);

			ID(const float x = sinf(_polar_angle) * sinf(_azimuth));
			ID(const float y = sinf(_polar_angle) * cosf(_azimuth));
			ID(const float z = cosf(_polar_angle));
			ID(const float xp =
			   sinf(p._polar_angle) * sinf(p._azimuth));
			ID(const float yp =
			   sinf(p._polar_angle) * cosf(p._azimuth));
			ID(const float zp = cosf(p._polar_angle));

			const float hav_dl =
				havf(p._polar_angle - _polar_angle);
			const float hav_dp =
				havf(p._azimuth - _azimuth);
			const float polar_metric =
				sinf(_polar_angle) * sinf(p._polar_angle);
			const float retval = 2.0F *
				asinf(sqrtf(hav_dp + polar_metric * hav_dl));

			I(retval >= 0 && retval <= M_PI);

			I(FEQ(retval, acosf(x * xp + y * yp + z * zp)));

			return (float)retval;
		}
    };

}

#endif // JETBASE_POLAR_H_
