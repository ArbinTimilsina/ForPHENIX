#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H
#include <jetbase/dbc.h>
#include <jetbase/specfunc.h>
#include <jetevent/lorentz.h>

namespace jet {

	/////////////////////////////////////////////////////////////////

	float abstract_lorentz_vector_t::
	max_pseudorapidity(const float perp, const float z)
	{
		static const float scale =
			1.0F - jet::machine_limit_t::spacing(perp);

		if(perp < 1.0F)
			return (z > 0 ? scale : -scale) *
				jet::machine_limit_t::max_pseudorapidity(perp);
		return (z > 0 ? scale : -scale) *
			asinhf(jet::machine_limit_t::max_magnitude(perp) *
				   scale / std::max(1.0F, perp));
	}

	float abstract_lorentz_vector_t::
	pseudorapidity_limit(const float __pseudorapidity,
						 const float perp)
	{
		// Rectangular masks for common ranges
		if((fabsf(__pseudorapidity) <= 8.0F &&
			perp < 2.2830403254864205452e+35F) ||
		   (fabsf(__pseudorapidity) <= 64.0F &&
			perp < 1.0914970017423396922e+11F))
		   return __pseudorapidity;

		static const float scale =
			1.0F - jet::machine_limit_t::spacing(perp);

		if(perp < 1.0F)
			return (__pseudorapidity > 0 ? scale : -scale) *
				jet::machine_limit_t::max_pseudorapidity(perp);

		const float limit = scale *
			asinhf(jet::machine_limit_t::max_magnitude(perp) *
				   scale / std::max(1.0F, perp));

		return __pseudorapidity > 0 ?
			std::min(limit, __pseudorapidity) :
			std::max(-limit, __pseudorapidity);
	}

	float abstract_lorentz_vector_t::perp_square(void) const
	{
		const float ret = x() * x() + y() * y();

		I(ret >= 0);

		return ret;
	}

	float abstract_lorentz_vector_t::perp(void) const
	{
		const float s = perp_square();
		float ret;

		I(s >= 0);

		ret = sqrtf(s);

		I(F(ret));
		I(ret >= fabsf(x()));
		I(ret >= fabsf(y()));
		I(ret <= fabsf(x()) + fabsf(y()));

		return ret;
	}

	float abstract_lorentz_vector_t::
	cartesian_magnitude_square(void) const
	{
		I(F(x()));
		I(F(y()));
		I(F(z()));

		const float s = x() * x() + y() * y() + z() * z();

		I(s >= 0);

		return s;
	}

	float abstract_lorentz_vector_t::cartesian_magnitude(void) const
	{
		const float s = cartesian_magnitude_square();

		I(s >= 0);

		const float m = sqrtf(s);

		I(m >= 0);
		I(FEQ(m * m, s));

		return m;
	}

	float abstract_lorentz_vector_t::mass_square(void) const
	{
		const float s = cartesian_magnitude_square();

		I(s >= 0);

		const float ret = time() * time() - s;

		// I(ret >= 0 || FEQ(ret, 0));

		return ret;
	}

	float abstract_lorentz_vector_t::mass(void) const
	{
		const float s = mass_square();

		// I(s >= 0 || FEQ(s, 0));

		// s can be slightly less than 0 due to cancellation.
		const float m = sqrtf(std::max(0.0F, s));

		I(m >= 0);

		return m;
	}

	float abstract_lorentz_vector_t::mass_perp_square(void) const
	{
		const float s = mass_square();

		// I(s >= 0 || FEQ(s, 0));

		const float t = perp_square();

		I(t >= 0 || FEQ(t, 0));

		return s + t;
	}

	float abstract_lorentz_vector_t::mass_perp(void) const
	{
		const float s = mass_perp_square();

		// s can be slightly less than 0 due to cancellation.
		const float m = sqrtf(std::max(0.0F, s));

		I(m >= 0);

		return m;
	}

	float abstract_lorentz_vector_t::pseudorapidity(void) const
	{
		const float s = perp_square();
		const float magnitude_square = s + z() * z();
		float ret = 0;

		I(magnitude_square >= 0);

		if(std::fpclassify(s) == FP_ZERO)
			ret = max_pseudorapidity(ret, z());
		else {		// std::fpclassify(s) != FP_ZERO
			// The accuracy crossover between the two algorithms in
			// the IEEE 754 single precision is dominated by the step
			// function of ulp(1 - s^2) ~ 2^(-24 - theta(x -
			// 1/sqrt(2))) around sin_polar_angle_square := s^2 =
			// 1/sqrt(2), which corresponds to eta == log(3 + 2 *
			// sqrt(2)) / 2 == 0.881373587.
			if(s < (float)M_SQRT1_2 * magnitude_square) {
				I(std::fpclassify(z()) != FP_ZERO);

				const float half_sign = z() > 0 ? 0.5F : -0.5F;
				const float sin_polar_angle_square =
					s / magnitude_square;

				// The expression atanh(sqrt(1 - s^2)) can be
				// rewritten into a form that is analytic for 0 <= s <
				// 1 by using the identity atanh(x) = (1 / 2) * log((1
				// + x) / (1 - x)), thus separating its log(2 / s)
				// pole around st == 0.
				ret = half_sign *
					logf((2.0F - sin_polar_angle_square +
						  2.0F * sqrtf(1.0F -
									   sin_polar_angle_square)) /
						 sin_polar_angle_square);
			}
			else {
				const float cos_polar_angle =
					std::fpclassify(magnitude_square) == FP_ZERO ?
					1.0F : z() / sqrtf(magnitude_square);

				I(FRANGE(cos_polar_angle, -1, 1));

				ret = atanhf(cos_polar_angle);

				I(FEQ(tanh(ret), cos_polar_angle));
				I(F(ret));
			}
			ret = pseudorapidity_limit(ret, perp());
		}

		return ret;
	}

	float abstract_lorentz_vector_t::rapidity(void) const
	{
		// This is likely not numerically stable

		return 0.5F * logf((time() + z()) / (time() - z()));
	}

	float abstract_lorentz_vector_t::sin_polar_angle(void) const
	{
		// Is this implementation the numerically best?
		const float t = perp() / z();

		I(t >= 0);

		const float ret = 1.0F / sqrtf(std::max(0.0F, 1.0F + t * t));

		I(FRANGE(ret, 0, 1));

		return ret;
	}

	float abstract_lorentz_vector_t::time_perp(void) const
	{
		const float s = sin_polar_angle();
		const float t = time();

		I(FRANGE(t, 0, 1));

		const float ret = s * t;

		I(ret <= t);

		return ret;
	}

	float abstract_lorentz_vector_t::azimuth(void) const
	{
		const float ret = (std::fpclassify(x()) == FP_ZERO &&
						   std::fpclassify(y()) == FP_ZERO) ?
			0.0F : atan2f(y(), x());

		I(F(ret));

		return ret;
	}

	float abstract_lorentz_vector_t::
	cartesian_dot(const abstract_lorentz_vector_t &v) const
	{
		const float s = x() * v.x() + y() * v.y() + z() * v.z();

		return s;
	}

	float abstract_lorentz_vector_t::
	pseudorapidity_edge_distance(
		const float pseudorapidity_min,
		const float pseudorapidity_max) const
	{
		const float eta = pseudorapidity();

		return std::min(
			eta - pseudorapidity_min,
			pseudorapidity_max - eta);
	}

	float abstract_lorentz_vector_t::
	azimuth_edge_distance(
		const float azimuth_min,
		const float azimuth_max) const
	{
		const float phi = azimuth();

		return std::min(
			jet::angular_range_reduce(phi - azimuth_min),
			jet::angular_range_reduce(azimuth_max - phi));
	}

	float abstract_lorentz_vector_t::
	azimuth_edge_distance(
		const float azimuth_min_0,
		const float azimuth_max_0,
		const float azimuth_min_1,
		const float azimuth_max_1) const
	{
		const float phi = azimuth();

		return std::max(
			std::min(
				jet::angular_range_reduce(phi - azimuth_min_0),
				jet::angular_range_reduce(azimuth_max_0 - phi)),
			std::min(
				jet::angular_range_reduce(phi - azimuth_min_1),
				jet::angular_range_reduce(azimuth_max_1 - phi)));
	}

	/////////////////////////////////////////////////////////////////

	lorentz_vector_t::
	lorentz_vector_t(const float time, const float x, const float y,
					 const float z)
	{
		// Testing for I(F({time, x, y, z})) is not possible here,
		// NAN is used for unknown vertices
		_x[0] = std::min(
			std::max(time, -machine_limit_t::max_cartesian(time)),
			machine_limit_t::max_cartesian(time));
		_x[1] = std::min(
			std::max(x, -machine_limit_t::max_cartesian(x)),
			machine_limit_t::max_cartesian(x));
		_x[2] = std::min(
			std::max(y, -machine_limit_t::max_cartesian(y)),
			machine_limit_t::max_cartesian(y));
		_x[3] = std::min(
			std::max(z, -machine_limit_t::max_cartesian(z)),
			machine_limit_t::max_cartesian(z));
	}

	lorentz_vector_t::lorentz_vector_t(const float time,
									   const float x[3])
	{
		// Testing for I(F({time, x, y, z})) is not possible here,
		// NAN is used for unknown vertices
		_x[0] = std::min(
			std::max(time, -machine_limit_t::max_cartesian(time)),
			machine_limit_t::max_cartesian(time));

		I(x != NULL);

		_x[1] = std::min(
			std::max(x[0], -machine_limit_t::max_cartesian(x[0])),
			machine_limit_t::max_cartesian(x[0]));
		_x[2] = std::min(
			std::max(x[1], -machine_limit_t::max_cartesian(x[1])),
			machine_limit_t::max_cartesian(x[1]));
		_x[3] = std::min(
			std::max(x[2], -machine_limit_t::max_cartesian(x[2])),
			machine_limit_t::max_cartesian(x[2]));
	}

	lorentz_vector_t lorentz_vector_t::
	transverse(const lorentz_vector_t &v) const
	{
		const float s = cartesian_dot(v) /
			v.cartesian_magnitude_square();

		I(F(s));

		const lorentz_vector_t ret = *this - v * s;

		return ret;
	}

	float lorentz_vector_t::
	transverse_magnitude(const lorentz_vector_t &v) const
	{
		const lorentz_vector_t t = transverse(v);

		I(FEQ(0.01 * t.cartesian_dot(v), 0));

		return t.cartesian_magnitude();
	}

}
