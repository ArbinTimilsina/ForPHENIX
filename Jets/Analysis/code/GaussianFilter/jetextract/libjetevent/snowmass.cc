#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H
#include <jetevent/snowmass.h>

// JET RECONSTRUCTION NUMERICS - SNOWMASS KINEMATICS

namespace jet {

	/////////////////////////////////////////////////////////////////

	snowmass_vector_t::
	snowmass_vector_t(const float energy, const float perp,
					  const float pseudorapidity,
					  const float azimuth)
	{
		// No precondition needed for energy
		_x[0] = energy;

		// Do not specify I(perp >= 0), as we must allow calorimeter
		// clusters with undefined transverse momentum

		_x[1] = perp;

		// I(F(pseudorapidity));

		_x[2] = pseudorapidity;
		_x[3] = angular_range_reduce(azimuth);

		I(FRANGE(_x[3], -M_PI, M_PI));
	}

	snowmass_vector_t::operator lorentz_vector_t(void) const
	{
		const float __time = time();
		const float __x = x();
		const float __y = y();
		const float __z = z();

		// I(F(__time));
		I(F(__x));
		I(F(__y));
		I(F(__z));

		return lorentz_vector_t(__time, __x, __y, __z);
	}

	void snowmass_vector_t::set_time(const float __time)
	{
		I(_x != NULL);

		// Do not specify I(energy >= 0), as we must allow tracks with
		// undefined energy

		_x[0] = __time;

		IG(FEQ(time(), __time),
		   !(std::isnan(time()) && std::isnan(__time)));
	}

	void snowmass_vector_t::set_perp_pseudorapidity_azimuth(
		const float perp, const float pseudorapidity,
		const float azimuth)
	{
		I(_x != NULL);
		I(F(perp));
		I(F(pseudorapidity));
		I(F(azimuth));
		I(FRANGE(perp, -machine_limit_t::max_cartesian(perp),
				 machine_limit_t::max_cartesian(perp)));

		_x[1] = perp;
		_x[2] = pseudorapidity;
		_x[3] = angular_range_reduce(azimuth);

		I(FEQ(sinf(_x[3]), sinf(azimuth)));
		I(FEQ(cosf(_x[3]), cosf(azimuth)));
		I(FRANGE(_x[3], -M_PI, M_PI));
	}

	void snowmass_vector_t::
	set_x_y_z(const float x, const float y, const float z)
	{
		I(F(x));
		I(F(y));
		I(F(z));
		I(FRANGE(x, -machine_limit_t::max_cartesian(x),
				 machine_limit_t::max_cartesian(x)));
		I(FRANGE(y, -machine_limit_t::max_cartesian(y),
				 machine_limit_t::max_cartesian(y)));
		I(FRANGE(z, -machine_limit_t::max_cartesian(z),
				 machine_limit_t::max_cartesian(z)));
#if defined(HAVE_ROOT) && !defined(NVERIFY)
		ID(TLorentzVector v = TLorentzVector(x, y, z, 0));
#endif // defined(HAVE_ROOT) && !defined(NVERIFY)

		// Transverse momentum

		const float perp_square = x * x + y * y;

		I(perp_square >= 0);

		_x[1] = sqrtf(perp_square);

		// I(F(_x[1]));
		I(_x[1] >= fabs(x));
		I(_x[1] >= fabs(y));
		I(_x[1] <= fabs(x) + fabs(y));

		// Pseudorapidity

		const float magnitude_square = perp_square + z * z;

		I(magnitude_square >= 0);

		if(std::fpclassify(perp_square) == FP_ZERO)
			_x[2] = max_pseudorapidity(_x[1], z);
		else {		// std::fpclassify(perp_square) != FP_ZERO
			// The accuracy crossover between the two algorithms in
			// the IEEE 754 single precision is dominated by the step
			// function of ulp(1 - s^2) ~ 2^(-24 - theta(x -
			// 1/sqrt(2))) around sin_polar_angle_square := s^2 =
			// 1/sqrt(2), which corresponds to eta == log(3 + 2 *
			// sqrt(2)) / 2 == 0.881373587.
			if(perp_square < (float)M_SQRT1_2 * magnitude_square) {
				I(std::fpclassify(z) != FP_ZERO);

				const float half_sign = z > 0 ? 0.5F : -0.5F;
				const float sin_polar_angle_square =
					perp_square / magnitude_square;

				// The expression atanh(sqrt(1 - s^2)) can be
				// rewritten into a form that is analytic for 0 <= s <
				// 1 by using the identity atanh(x) = (1 / 2) * log((1
				// + x) / (1 - x)), thus separating its log(2 / s)
				// pole around st == 0.
				_x[2] = half_sign *
					logf((2.0F - sin_polar_angle_square +
						  2.0F * sqrtf(1.0F -
									   sin_polar_angle_square)) /
						 sin_polar_angle_square);
			}
			else {
				const float cos_polar_angle =
					std::fpclassify(magnitude_square) == FP_ZERO ?
					1.0F : z / sqrtf(magnitude_square);

				I(FRANGE(cos_polar_angle, -1, 1));

				_x[2] = atanhf(cos_polar_angle);

				I(FEQ(tanh(_x[2]), cos_polar_angle));
				I(F(_x[2]));
			}
			_x[2] = pseudorapidity_limit(_x[2], _x[1]);
		}

		// Azimuth

		_x[3] = (std::fpclassify(x) == FP_ZERO &&
				 std::fpclassify(y) == FP_ZERO) ?
			0.0F : atan2f(y, x);

		I(F(_x[3]));
#if defined(HAVE_ROOT) && !defined(NVERIFY)
		// No postcondition needed for energy (since there was no
		// precondition)
		I(FEQ(_x[1], v.Pt()));
		IG(FEQ(_x[2], v.Eta()), _x[2] >= -8 && _x[2] <= 8);
		I(FEQ(_x[3], v.Phi()));
#endif // defined(HAVE_ROOT) && !defined(NVERIFY)
		I(FEQ(this->x(), x));
		I(FEQ(this->y(), y));
	}

	float snowmass_vector_t::x(void) const
	{
		I(_x != NULL);
#if defined(HAVE_ROOT) && !defined(NVERIFY)
		ID(TLorentzVector v);
		IS(v.SetPtEtaPhiE(_x[1], _x[2], _x[3], _x[0]));
#endif // defined(HAVE_ROOT) && !defined(NVERIFY)

		const float __perp = perp();
		const float __azimuth = azimuth();

		//I(__perp >= 0);

		const float ret = __perp * cosf(__azimuth);

#if defined(HAVE_ROOT) && !defined(NVERIFY)
		IG(FEQ(ret, v.Vect().Px()), _x[2] >= -8 && _x[2] <= 8 && __perp >= 0);
#endif // defined(HAVE_ROOT) && !defined(NVERIFY)
		I(F(ret));
		I(FRANGE(ret, -machine_limit_t::max_cartesian(ret),
				 machine_limit_t::max_cartesian(ret)));

		return ret;
	}

	float snowmass_vector_t::y(void) const
	{
		I(_x != NULL);
#if defined(HAVE_ROOT) && !defined(NVERIFY)
		ID(TLorentzVector v);
		IS(v.SetPtEtaPhiE(_x[1], _x[2], _x[3], _x[0]));
#endif // defined(HAVE_ROOT) && !defined(NVERIFY)

		const float __perp = perp();
		const float __azimuth = azimuth();

		//I(__perp >= 0);

		const float ret = __perp * sinf(__azimuth);

#if defined(HAVE_ROOT) && !defined(NVERIFY)
		IG(FEQ(ret, v.Vect().Py()), _x[2] >= -8 && _x[2] <= 8 && __perp >= 0);
#endif // defined(HAVE_ROOT) && !defined(NVERIFY)
		I(F(ret));
		I(FRANGE(ret, -machine_limit_t::max_cartesian(ret),
				 machine_limit_t::max_cartesian(ret)));

		return ret;
	}

	float snowmass_vector_t::z(void) const
	{
		I(_x != NULL);
#if defined(HAVE_ROOT) && !defined(NVERIFY)
		ID(TLorentzVector v);
		IS(v.SetPtEtaPhiE(_x[1], _x[2], _x[3], _x[0]));
#endif // defined(HAVE_ROOT) && !defined(NVERIFY)

		const float __perp = perp();
		const float __pseudorapidity = pseudorapidity();

		//I(__perp >= 0);

		// Note: cot(2 * atan(exp(-t))) == sinh(t)
		const float ret = __perp * sinhf(__pseudorapidity);

#if 0
		if(!F(ret)) {
			fprintf(stderr, "__perp = %.8g, __pseudorapidity = %.8g, "
					"sinhf(__pseudorapidity) = %.8g\n",
					__perp, __pseudorapidity,
					sinhf(__pseudorapidity));
		}
#endif
		IG(F(ret), __perp >= 0.0F);

		return ret;
	}

	// Calculate the Cartesian norm (|p|)
	float snowmass_vector_t::cartesian_magnitude(void) const
	{
		I(_x != NULL);
		// I(F(_x[1]));
		I(F(_x[2]));
		I(F(_x[3]));
		// Note: IG(_x[1] >= 0, std::isfinite(_x[1])) cannot
		// necessarily be satisfied with heavy ion events
#if defined(HAVE_ROOT) && !defined(NVERIFY)
		ID(TLorentzVector v);
		IS(v.SetPtEtaPhiE(_x[1], _x[2], _x[3], _x[0]));
#endif // defined(HAVE_ROOT) && !defined(NVERIFY)

		const float perp = _x[1], eta = _x[2];
		// Note: sin(2 * atan(exp(-t))) == sech(t)
		const float s = coshf(eta);

		I(s >= 0);

		const float ret = perp * s;

#if defined(HAVE_ROOT) && !defined(NVERIFY)
		IG(FEQ(ret, v.Vect().Mag()), _x[1] >= 0 && _x[2] >= -8 && _x[2] <= 8);
#endif // defined(HAVE_ROOT) && !defined(NVERIFY)
		IG(ret >= 0, _x[1] >= 0);

		return ret;
	}

	// Calculate the Cartesian norm squared (|p|^2)
	float snowmass_vector_t::cartesian_magnitude_square(void) const
	{
		I(_x != NULL);
		// I(F(_x[1]));
		I(F(_x[2]));
		I(F(_x[3]));
#if defined(HAVE_ROOT) && !defined(NVERIFY)
		ID(TLorentzVector v);
		IS(v.SetPtEtaPhiE(_x[1], _x[2], _x[3], _x[0]));
#endif // defined(HAVE_ROOT) && !defined(NVERIFY)

		const float s = cartesian_magnitude();
		const float ret = s * s;

		// Note: ret can overflow even if s is finite
#if defined(HAVE_ROOT) && !defined(NVERIFY)
		IG(FEQ(ret, v.Vect().Mag2()), _x[2] >= -8 && _x[2] <= 8);
#endif // defined(HAVE_ROOT) && !defined(NVERIFY)
		IG(ret >= 0, F(_x[1]) && F(_x[2]) && F(_x[3]));

		return ret;
	}

	float snowmass_vector_t::cos_polar_angle(void) const
	{
		I(_x != NULL);
		I(F(_x[2]));
#if defined(HAVE_ROOT) && !defined(NVERIFY)
		ID(TLorentzVector v);
		IS(v.SetPtEtaPhiE(1, _x[2], _x[3], 0));
#endif // defined(HAVE_ROOT) && !defined(NVERIFY)

		const float ret = tanhf(_x[2]);

#if defined(HAVE_ROOT) && !defined(NVERIFY)
		I(FEQ(ret, v.CosTheta()));
#endif // defined(HAVE_ROOT) && !defined(NVERIFY)
		I(F(ret));
		I(FRANGE(ret, -1, 1));
		// Assertion: sign(ret) == sign(_x[2])
		I((ret > 0 ? 1 : ret < 0 ? -1 : 0) ==
		  (_x[2] > 0 ? 1 : _x[2] < 0 ? -1 : 0));

		return ret;
	}

	float snowmass_vector_t::sin_polar_angle(void) const
	{
		I(_x != NULL);
		I(F(_x[2]));

		const float ret = 1 / coshf(_x[2]);

		I(F(ret));
		I(FRANGE(ret, 0, 1));

		return ret;
	}

	// Arithmetics

	snowmass_vector_t
	snowmass_vector_t::recombine(const snowmass_vector_t &s) const
	{
		// Cody and Waite expansion of 2 pi
		static const float two_pi_0 =  6.283185482F;
		static const float two_pi_1 = -1.7484556001e-7F;

		// Transverse sum

		const float perp = _x[1] + s._x[1];

		I(perp > 0);

		// Pseudorapidity sum

		const float inv_perp = 1 / perp;
		const float eta =
			(_x[1] * _x[2] + s._x[1] * s._x[2]) * inv_perp;

		I(_x[3] > (float)(-M_PI) && _x[3] <= (float)M_PI);

		// Azimuthal sum

		const float dphi = _x[3] - s._x[3];
		float phi;

		// Handling of cyclic redundancy
		I(dphi >= -(float)(2 * M_PI) && dphi <= (float)(2 * M_PI));

		if(dphi < -M_PI)
			phi = (_x[1] * ((_x[3] + two_pi_0) + two_pi_1) +
				   s._x[1] * s._x[3]) * inv_perp;
		else if(dphi > M_PI)
			phi = (_x[1] * ((_x[3] - two_pi_0) - two_pi_1) +
				   s._x[1] * s._x[3]) * inv_perp;
		else
			phi = (_x[1] * _x[3] +
				   s._x[1] * s._x[3]) * inv_perp;

		return snowmass_vector_t(_x[0] + s._x[0], perp, eta, phi);
	}

	float snowmass_vector_t::
	cartesian_dot(const snowmass_vector_t &s) const
	{
		I(_x != NULL);
		// I(F(_x[1]));
		I(F(_x[2]));
		I(F(_x[3]));
		I(FRANGE(_x[3], -M_PI, M_PI));
		I(s._x != NULL);
		// I(F(s._x[1]));
		I(F(s._x[2]));
		I(F(s._x[3]));
		I(FRANGE(s._x[3], -M_PI, M_PI));

#if defined(HAVE_ROOT) && !defined(NVERIFY)
		TLorentzVector tlv1;
		TLorentzVector tlv2;
		tlv1.SetPtEtaPhiE(perp(), pseudorapidity(), azimuth(),
						  time());
		tlv2.SetPtEtaPhiE(s.perp(), s.pseudorapidity(), s.azimuth(),
						  s.time());
#endif // defined(HAVE_ROOT) && !defined(NVERIFY)

		// Note: cot(2 * atan(exp(-t))) == sinh(t)
		const float ret = _x[1] * s._x[1] *
			(cosf(_x[3] - s._x[3]) + sinhf(_x[2]) * sinhf(s._x[2]));

#if defined(HAVE_ROOT) && !defined(NVERIFY)
		IG(FEQ(ret, tlv1.Vect().Dot(tlv2.Vect())),
		   _x[1] >= 0 && _x[2] >= -2 && _x[2] <= 2 &&
		   s._x[1] >= 0 && s._x[2] >= -2 && s._x[2] <= 2);
#endif // defined(HAVE_ROOT) && !defined(NVERIFY)
		IG(F(ret),
		   F(_x[1]) && F(_x[2]) && F(_x[3]) &&
		   F(s._x[1]) && F(s._x[2]) && F(s._x[3]) &&
		   _x[2] >= -8 && _x[2] <= 8 &&
		   s._x[2] >= -8 && s._x[2] <= 8);
		ID(const float ret_bound =
		   cartesian_magnitude() * s.cartesian_magnitude());
		IG(FRANGE(ret, -ret_bound, ret_bound), 
		   F(_x[1]) && F(_x[2]) && F(_x[3]) &&
		   F(s._x[1]) && F(s._x[2]) && F(s._x[3]) &&
		   _x[1] >= 0 && _x[2] >= -8 && _x[2] <= 8 && 
		   s._x[1] >= 0 && s._x[2] >= -8 && s._x[2] <= 8);

		return ret;
	}

	float snowmass_vector_t::
	radial_distance_square(const snowmass_vector_t &v) const
	{
		const float pseudorapidity_difference =
			v.pseudorapidity() - pseudorapidity();
		const float azimuth_difference =
			angular_range_reduce(v.azimuth() - azimuth());

		I(F(pseudorapidity_difference));
		I(F(azimuth_difference));

		const float ret =
			pseudorapidity_difference *
			pseudorapidity_difference +
			azimuth_difference * azimuth_difference;

		I(ret >= 0);

		return ret;
	}

	float snowmass_vector_t::
	radial_distance(const snowmass_vector_t &v) const
	{
		const float square = radial_distance_square(v);

		I(square >= 0);

		const float ret = sqrtf(square);

		I(ret >= 0);

		return ret;
	}

	float snowmass_vector_t::radial_origin_distance(void) const
	{
		const float square = pseudorapidity() * pseudorapidity() +
			azimuth() * azimuth();

		I(square >= 0);

		const float ret = sqrtf(square);

		I(ret >= 0);

		return ret;
	}

	// Specialized routines for realistic detectors

	void snowmass_vector_t::set_lightlike_perp(void)
	{
		const float ct = cos_polar_angle();
		const float st = sqrtf(1 - ct * ct);

		_x[1] = _x[0] * st;

		// I(F(_x[1]));
		I(_x[1] <= _x[0]);
	}

	void snowmass_vector_t::set_lightlike_time(void)
	{
		const float ct = cos_polar_angle();
		const float st = sqrtf(1 - ct * ct);

		I(st >= 0);

		if(std::fpclassify(st) == FP_ZERO)
			_x[0] = machine_limit_t::max_cartesian(_x[0]);
		else
			_x[0] = _x[1] / st;

		I(F(_x[0]));
		//I(_x[0] >= _x[1]);
	}
}
