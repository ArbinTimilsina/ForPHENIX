#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif // HAVE_CONFIG_H

#include <cstdlib>
#include <jetrec/rec.h>

namespace jet {

	const char *reconstruction_filtering_t::status_str(
		const unsigned int status)
	{
		switch(status) {
		case STATUS_CONTINUE:
			return "STATUS_CONTINUE";
			break;
		case STATUS_CONVERGENCE:
			return "STATUS_CONVERGENCE";
			break;
		case STATUS_HESSIAN_INDEFINITE:
			return "STATUS_HESSIAN_INDEFINITE";
			break;
		case STATUS_HESSIAN_SINGULAR:
			return "STATUS_HESSIAN_SINGULAR";
			break;
		case STATUS_SADDLE_POINT:
			return "STATUS_SADDLE_POINT";
			break;
		case STATUS_NOT_NUMERICAL:
			return "STATUS_NOT_NUMERICAL";
			break;
		case STATUS_LINE_SEARCH_CONTINUE:
			return "STATUS_LINE_SEARCH_CONTINUE";
			break;
		case STATUS_LINE_SEARCH_CONVERGENCE:
			return "STATUS_LINE_SEARCH_CONVERGENCE";
			break;
		case STATUS_LINE_SEARCH_STEP_LIMIT:
			return "STATUS_LINE_SEARCH_STEP_LIMIT";
			break;
		case STATUS_LINE_SEARCH_INVALID_BRACKET:
			return "STATUS_LINE_SEARCH_INVALID_BRACKET";
			break;
		case STATUS_LINE_SEARCH_NOT_DESCENDING:
			return "STATUS_LINE_SEARCH_NOT_DESCENDING";
			break;
		case STATUS_LINE_SEARCH_INVALID_RANGE:
			return "STATUS_LINE_SEARCH_INVALID_RANGE";
			break;
		case STATUS_LINE_SEARCH_NOT_NUMERICAL:
			return "STATUS_LINE_SEARCH_NOT_NUMERICAL";
			break;
		default:
			return "(unknown status)";
		}
	}

	// The dimensions are 1: pseudorapidity, 2: azimuth
	// The gradient is stored as { g(1), g(2) }
	// The Hessian is stored as { h(1, 1), h(2, 2), h(1, 2) }

	void reconstruction_filtering_t::
	modify_hessian_spectral(
		float modified_hessian[], const float hessian[]) const
	{
		// Calculate the eigensystem of the Hessian
		const float d = hessian[0] - hessian[1];
		const float discriminant =
			d * d + 4.0F * hessian[2] * hessian[2];
		const float sqrt_discriminant = sqrtf(discriminant);
		const float s = hessian[0] + hessian[1];
		const float eigenvalue_1 = 0.5F * (s - sqrt_discriminant);
		const float eigenvalue_2 = 0.5F * (s + sqrt_discriminant);
		const float e1 =
			(d - sqrt_discriminant) / (2.0F * hessian[2]);
		const float e2 =
			(d + sqrt_discriminant) / (2.0F * hessian[2]);
		const float eigenvector_11 = e1 / sqrtf(e1 * e1 + 1.0F);
		const float eigenvector_12 = 1.0F / sqrtf(e1 * e1 + 1.0F);
		const float eigenvector_21 = e2 / sqrtf(e2 * e2 + 1.0F);
		const float eigenvector_22 = 1.0F / sqrtf(e2 * e2 + 1.0F);
		// Enforce a negative definite Hessian by flipping the
		// positive eigenvalues
		const float modified_eigenvalue_1 = -fabsf(eigenvalue_1);
		const float modified_eigenvalue_2 = -fabsf(eigenvalue_2);
		// Unitary transformation to obtain the modified Hessian

		modified_hessian[0] =
			modified_eigenvalue_1 *
			eigenvector_11 * eigenvector_11 +
			modified_eigenvalue_2 *
			eigenvector_12 * eigenvector_12;
		modified_hessian[2] =
			modified_eigenvalue_1 *
			eigenvector_11 * eigenvector_21 +
			modified_eigenvalue_2 *
			eigenvector_12 * eigenvector_22;
		modified_hessian[1] =
			modified_eigenvalue_1 *
			eigenvector_21 * eigenvector_21 +
			modified_eigenvalue_2 *
			eigenvector_22 * eigenvector_22;
	}

	void reconstruction_filtering_t::
	modify_hessian_modified_cholesky(
		float modified_hessian[], const float hessian[]) const
	{
		static const float inverse_nu = 0.57735027F; // 1 / sqrt(3)
		static const float epsilon = sqrtf(FLT_EPSILON);
		const float abs_hessian_0 = fabsf(hessian[0]);
		const float abs_hessian_1 = fabsf(hessian[1]);
		const float abs_hessian_2 = fabsf(hessian[2]);
		const float beta_square =
			std::max(std::max(abs_hessian_0, abs_hessian_1),
					 abs_hessian_2 * inverse_nu);
		const float hessian_2_square = hessian[2] * hessian[2];

		// Modified Cholesky factorization of the hessian into LDL^T
		// (elements denoted as * are not stored)
		//
		// L = { { 1, 0 }, { *, 1 } }
		// D = diag({ modified_hessian[0], * })
		// E = diag({ modified_hessian[0] - hessian[0],
		// modified_hessian[1] - hessian[1] })

		modified_hessian[0] =
			-std::max(std::max(epsilon, abs_hessian_0),
					  hessian_2_square / beta_square);
		modified_hessian[1] =
			hessian[1] - hessian_2_square / modified_hessian[0];
		modified_hessian[1] = hessian[1] -
			std::max(epsilon, fabsf(modified_hessian[1])) -
			modified_hessian[1];
		modified_hessian[2] = hessian[2];
	}

	void reconstruction_filtering_t::
	constrain_newton_direction_pseudorapidity(
		float &newton_direction_pseudorapidity,
		float &newton_direction_azimuth,
		const float pseudorapidity) const
	{
		const float result_pseudorapidity =
			pseudorapidity + newton_direction_pseudorapidity;

		newton_direction_azimuth = angular_range_reduce(
			newton_direction_azimuth);
		if(result_pseudorapidity < _pseudorapidity_range.first) {
			const float new_direction_pseudorapidity =
				_pseudorapidity_range.first - pseudorapidity;

			newton_direction_azimuth =
				newton_direction_azimuth *
				new_direction_pseudorapidity /
				newton_direction_pseudorapidity;
			newton_direction_pseudorapidity =
				new_direction_pseudorapidity;
		}
		else if(result_pseudorapidity > _pseudorapidity_range.second) {
			const float new_direction_pseudorapidity =
				_pseudorapidity_range.second - pseudorapidity;

			newton_direction_azimuth =
				newton_direction_azimuth *
				new_direction_pseudorapidity /
				newton_direction_pseudorapidity;
			newton_direction_pseudorapidity =
				new_direction_pseudorapidity;
		}
	}

	float reconstruction_filtering_t::evaluate_quadratic(
		const float x0, const float x1,
		const float gradient[], const float hessian[]) const
	{
		// Quadratic model f(x + s) - f(x)
		return gradient[0] * x0 + gradient[1] * x1 + 0.5F *
			(hessian[0] * x0 * x0 + x1 *
			 (2.0F * hessian[2] * x0 + hessian[1] * x1));
	}

	float reconstruction_filtering_t::determinant_2(const float a[])
		const
	{
		return a[1] * a[0] - a[2] * a[2];
	}

	void reconstruction_filtering_t::solve_quadratic_unconstrained_2(
		float &x0, float &x1,
		const float gradient[], const float hessian[],
		const float determinant) const
	{
		x0 = (hessian[2] * gradient[1] -
			  hessian[1] * gradient[0]) / determinant;
		x1 = (hessian[2] * gradient[0] -
			  hessian[0] * gradient[1]) / determinant;
	}

	bool reconstruction_filtering_t::solve_quadratic_constrained_2_0(
		float &x0, const float x1,
		const float x0_min, const float x0_max,
		const float gradient[], const float hessian[]) const
	{
		x0 = -(gradient[0] + hessian[2] * x1) / hessian[0];
		return x0 >= x0_min && x0 <= x0_max;
	}

	bool reconstruction_filtering_t::solve_quadratic_constrained_2_1(
		const float x0, float &x1,
		const float x1_min, const float x1_max,
		const float gradient[], const float hessian[]) const
	{
		x1 = -(gradient[1] + hessian[2] * x0) / hessian[1];
		return x1 >= x1_min && x1 <= x1_max;
	}

	void reconstruction_filtering_t::solve_quadratic_constrained_2(
		float &max_x0, float &max_x1,
		const float x0_min, const float x0_max,
		const float x1_min, const float x1_max,
		const float gradient[], const float hessian[]) const
	{
		// Active set is not empty
		float active_x0;
		float active_x1;
		bool no_active = true;
		float max_q = NAN;
		float q;

		// 1. x0, lower
		active_x0 = x0_min;
		if(solve_quadratic_constrained_2_1(
			active_x0, active_x1, x1_min, x1_max,
			gradient, hessian)) {
			q = evaluate_quadratic(
				active_x0, active_x1,
				gradient, hessian);

			if(no_active || q > max_q) {
				max_q = q;
				max_x0 = active_x0;
				max_x1 = active_x1;
				no_active = false;
			}
		}
		// 2. x1, lower
		active_x1 = x1_min;
		if(solve_quadratic_constrained_2_0(
			active_x0, active_x1, x0_min, x0_max,
			gradient, hessian)) {
			q = evaluate_quadratic(
				active_x0, active_x1,
				gradient, hessian);

			if(no_active || q > max_q) {
				max_q = q;
				max_x0 = active_x0;
				max_x1 = active_x1;
				no_active = false;
			}
		}
		// 3. x0, upper
		active_x0 = x0_max;
		if(solve_quadratic_constrained_2_1(
			active_x0, active_x1, x1_min, x1_max,
			gradient, hessian)) {
			q = evaluate_quadratic(
				active_x0, active_x1,
				gradient, hessian);

			if(no_active || q > max_q) {
				max_q = q;
				max_x0 = active_x0;
				max_x1 = active_x1;
				no_active = false;
			}
		}
		// 4. x1, upper
		active_x1 = x1_max;
		if(solve_quadratic_constrained_2_0(
			active_x0, active_x1, x0_min, x0_max,
			gradient, hessian)) {
			q = evaluate_quadratic(
				active_x0, active_x1,
				gradient, hessian);

			if(no_active || q > max_q) {
				max_q = q;
				max_x0 = active_x0;
				max_x1 = active_x1;
				no_active = false;
			}
		}
		// 5. Lower left
		active_x0 = x0_min;
		active_x1 = x1_min;
		q = evaluate_quadratic(
			active_x0, active_x1, gradient, hessian);
		if(no_active || q > max_q) {
			max_q = q;
			max_x0 = active_x0;
			max_x1 = active_x1;
			no_active = false;
		}
		// 6. Upper left
		active_x0 = x0_min;
		active_x1 = x1_max;
		q = evaluate_quadratic(
			active_x0, active_x1, gradient, hessian);
		if(no_active || q > max_q) {
			max_q = q;
			max_x0 = active_x0;
			max_x1 = active_x1;
			no_active = false;
		}
		// 7. Lower right
		active_x0 = x0_max;
		active_x1 = x1_min;
		q = evaluate_quadratic(
			active_x0, active_x1, gradient, hessian);
		if(no_active || q > max_q) {
			max_q = q;
			max_x0 = active_x0;
			max_x1 = active_x1;
			no_active = false;
		}
		// 8. Upper right
		active_x0 = x0_max;
		active_x1 = x1_max;
		q = evaluate_quadratic(
			active_x0, active_x1, gradient, hessian);
		if(no_active || q > max_q) {
			max_q = q;
			max_x0 = active_x0;
			max_x1 = active_x1;
			no_active = false;
		}
#if 1
		if(!(max_q >= 0)) {
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": information: constrained solver failed for "
				"x = {"
					  << mathematica_form(max_x0) << ','
					  << mathematica_form(max_x1)
					  << ", gradient = {"
					  << mathematica_form(gradient[0]) << ','
					  << mathematica_form(gradient[1])
					  << "}, hessian = {"
					  << mathematica_form(hessian[0]) << ','
					  << mathematica_form(hessian[0]) << ','
					  << mathematica_form(hessian[0]) << '}'
					  << std::endl;
		}

		I(max_q >= 0);
		IG(max_q != 0,
		   x0_min > 0 && x0_max < 0 && x1_min > 0 && x1_max < 0 &&
		   gradient[0] != 0 && gradient[1] != 0 &&
		   hessian[0] != 0 && hessian[1] != 0 && hessian[2] != 0);
#endif
	}

	int reconstruction_filtering_t::
	line_search_more_thuente_step(
		float &step_best, float &function_best,
		float &derivative_best,
		float &step_second, float &function_second,
		float &derivative_second,
		float &step_current, float &function_current,
		float &derivative_current, bool &bracketed,
		const float step_min, const float step_max) const
	{
		// Note: Understanding the following code requires extensive
		// familiarity with line search algorithms and reading of
		// Mor\'e & Thuente (1994).

		// Check the input parameters for errors
		if(!(std::isfinite(function_current) &&
			 std::isfinite(function_best))) {
			return STATUS_LINE_SEARCH_NOT_NUMERICAL;
		}
		if((bracketed &&
			(step_current <= std::min(step_best, step_second) ||
			 step_current >= std::max(step_best, step_second)))) {
#if 0
			std::cerr << __FILE__ << ':' << __LINE__ << ": information: "
				"invalid bracketed step configuration "
				"(step_current = "
					  << mathematica_form(step_current)
					  << "; step_best = "
					  << mathematica_form(step_best)
					  << "; step_second = "
					  << mathematica_form(step_second) << ";)"
					  << std::endl;
#endif
			return STATUS_LINE_SEARCH_INVALID_BRACKET;
		}
		else if(derivative_best * (step_current - step_best) >= 0) {
			std::cerr << __FILE__ << ':' << __LINE__ << ": information: "
				"not a descending direction (derivative_best = "
					  << mathematica_form(derivative_best)
					  << "; step_best = "
					  << mathematica_form(step_best)
					  << "; step_current = "
					  << mathematica_form(step_current) << ";)"
					  << std::endl;
			return STATUS_LINE_SEARCH_NOT_DESCENDING;
		}
		else if(step_max < step_min) {
				std::cerr << __FILE__ << ':' << __LINE__
						  << ": information: invalid step range {"
						  << mathematica_form(step_min) << ','
						  << mathematica_form(step_max) << '}'
						  << std::endl;
			return STATUS_LINE_SEARCH_INVALID_RANGE;
		}

		// Determine if the derivatives have opposite signs
		const int signbit_derivative_current =
			std::signbit(derivative_current);
		const int signbit_derivative_best =
			std::signbit(derivative_best);
		const bool opposite_sign_derivatives =
			(signbit_derivative_current ^
			 signbit_derivative_best) != 0;
		bool bound;
		float step_trial;

		// Case 1 of Mor\'e & Thuente (1994) (on p. 299): f_t > f_l (a
		// higher function value).
		//
		// The minimum is bracketed. If the cubic step is closer to
		// step0 than the quadratic step, the cubic step is taken,
		// else the average of the cubic and quadratic steps is taken.
		if(function_current > function_best) {
			bound = true;

			const float step_current_step_best =
				step_current - step_best;
			const float theta =
				(function_best - function_current) * 3 /
				step_current_step_best +
				derivative_best + derivative_current;
			const float fabs_theta = fabsf(theta);
			const float fabs_derivative_best =
				fabsf(derivative_best);
			const float fabs_derivative_current =
				fabsf(derivative_current);
			const float s =
				std::max(std::max(fabs_theta,
								  fabs_derivative_best),
						 fabs_derivative_current);
			const float theta_s = theta / s;
			float gamma =
				s * sqrtf(theta_s * theta_s - (derivative_best / s) *
						  (derivative_current / s));

			if(step_current < step_best) {
				gamma = -gamma;
			}

			const float p = gamma - derivative_best + theta;
			const float q = gamma - derivative_best + gamma +
				derivative_current;
			const float r = p / q;
			const float step_cubic =
				step_best + r * step_current_step_best;
			const float step_quadratic = step_best +
				derivative_best /
				  ((function_best - function_current) /
				   step_current_step_best +
				   derivative_best) / 2 *
				step_current_step_best;
			const float fabs_step_cubic_step_best =
				fabsf(step_cubic - step_best);
			const float fabs_step_quadratic_step_best =
				fabsf(step_quadratic - step_best);

			step_trial = (fabs_step_cubic_step_best <
						  fabs_step_quadratic_step_best) ?
				step_cubic :
				step_cubic + (step_quadratic - step_cubic) / 2;
			bracketed = true;
		}
		// Case 2 of Mor\'e & Thuente (1994) (on p. 299): f_t <= f_l
		// and g_t * g_l < 0 (a lower function value and derivatives
		// of opposite sign).
		//
		// The minimum is bracketed. If the cubic step is closer to
		// step_best than the quadratic (secant) step, the cubic step
		// is taken, else the quadratic step is taken.
		else if(opposite_sign_derivatives) {
			bound = false;

			const float step_current_step_best =
				step_current - step_best;
			const float theta =
				3 * (function_best - function_current) /
				step_current_step_best +
				derivative_best + derivative_current;
			const float fabs_theta = fabsf(theta);
			const float fabs_derivative_best =
				fabsf(derivative_best);
			const float fabs_derivative_current =
				fabsf(derivative_current);
			const float s =
				std::max(std::max(fabs_theta,
								  fabs_derivative_best),
						 fabs_derivative_current);
			const float theta_s = theta / s;
			float gamma =
				s * sqrtf(theta_s * theta_s - (derivative_best / s) *
						  (derivative_current / s));

			if(step_current > step_best) {
				gamma = -gamma;
			}

			const float p = gamma - derivative_current + theta;
			const float q = gamma - derivative_current + gamma +
				derivative_best;
			const float r = p / q;

			const float step_cubic =
				step_current - r * step_current_step_best;
			const float step_quadratic = step_current -
				(derivative_current /
				 (derivative_current - derivative_best)) *
				step_current_step_best;
			const float fabs_step_cubic_step_current =
				fabsf(step_cubic - step_current);
			const float fabs_step_quadratic_step_current =
				fabsf(step_quadratic - step_current);

			step_trial = (fabs_step_cubic_step_current >
						  fabs_step_quadratic_step_current) ?
				step_cubic : step_quadratic;
			bracketed = true;
		}
		// Case 3 of Mor\'e & Thuente (1994) (on p. 299f): f_t <= f_l,
		// g_t * g_l >= 0, and |g_t| <= |g_l| (a lower function value,
		// derivatives of the same sign, and the magnitude of the
		// derivative decreases).
		//
		// The cubic step is only used if the cubic tends to infinity
		// in the direction of the step or if the minimum of the cubic
		// is beyond step. Otherwise the cubic step is defined to be
		// either step_min or step_max. The quadratic (secant) step is
		// also computed and if the minimum is bracketed then the the
		// step closest to step_best is taken, else the step farthest
		// away is taken.
		else if(fabsf(derivative_current) < fabsf(derivative_best)) {
			bound = true;

			const float step_current_step_best =
				step_current - step_best;
			const float theta =
				3 * (function_best - function_current) /
				step_current_step_best +
				derivative_best + derivative_current;
			const float fabs_theta = fabsf(theta);
			const float fabs_derivative_best =
				fabsf(derivative_best);
			const float fabs_derivative_current =
				fabsf(derivative_current);
			const float s =
				std::max(std::max(fabs_theta,
								  fabs_derivative_best),
						 fabs_derivative_current);
			// The case gamma = 0 only arises if the cubic does not
			// tend to infinity in the direction of the step.
			const float theta_s = theta / s;
			const float g = theta_s * theta_s -
				(derivative_best / s) * (derivative_current / s);
			const float zero = 0;
			float gamma =
				s * sqrtf(std::max(zero, g));

			if(step_current > step_best) {
				gamma = -gamma;
			}

			const float p = (gamma - derivative_current) + theta;
			const float q =
				(gamma + (derivative_best - derivative_current)) +
				gamma;
			const float r = p / q;
			float step_cubic;

			if(r < 0 && std::fpclassify(gamma) != FP_ZERO) {
				step_cubic =
					step_current - r * step_current_step_best;
			}
			else if(step_current > step_best) {
				step_cubic = step_max;
			}
			else {
				step_cubic = step_min;
			}

			const float step_quadratic = step_current -
				(derivative_current /
				 (derivative_current - derivative_best)) *
				step_current_step_best;
			const float fabs_step_current_step_cubic =
				fabsf(step_current - step_cubic);
			const float fabs_step_current_step_quadratic =
				fabsf(step_current - step_quadratic);

			step_trial = bracketed ?
				(fabs_step_current_step_cubic <
				 fabs_step_current_step_quadratic ?
				 step_cubic : step_quadratic) :
				(fabs_step_current_step_cubic >
				 fabs_step_current_step_quadratic ?
				 step_cubic : step_quadratic);
		}
		// Case 3 of Mor\'e & Thuente (1994) (on p. 300): f_t <= f_l,
		// g_t * g_l >= 0, and |g_t| > |g_l| (a lower function value,
		// derivatives of the same sign, and the magnitude of the
		// derivative does not decrease).
		//
		// If the minimum is not bracketed, the step is either
		// step_min or step_max, else the cubic step is taken.
		else {
			bound = false;
			if(bracketed) {
				const float step_second_step_current =
					step_second - step_current;
				const float theta =
					3 * (function_current - function_second) /
					step_second_step_current +
					derivative_second + derivative_current;
				const float fabs_theta = fabsf(theta);
				const float fabs_derivative_second =
					fabsf(derivative_second);
				const float fabs_derivative_current =
					fabsf(derivative_current);
				const float s =
					std::max(std::max(fabs_theta,
									  fabs_derivative_second),
							 fabs_derivative_current);
				const float theta_s = theta / s;
				float gamma =
					s * sqrtf(theta_s * theta_s -
							  (derivative_second / s) *
							  (derivative_current / s));

				if(step_current > step_second) {
					gamma = -gamma;
				}

				const float p = (gamma - derivative_current) +
					theta;
				const float q = ((gamma - derivative_current) +
								  gamma) + derivative_second;
				const float r = p / q;
				const float step_cubic =
					step_current + r * step_second_step_current;

				step_trial = step_cubic;
			}
			else if(step_current > step_best) {
				step_trial = step_max;
			}
			else {
				step_trial = step_min;
			}
		}

		// Mor\'e & Thuente (1994), p. 297, "Modified Updating
		// Algorithm"
		//
		// Update the interval of uncertainty. This update does not
		// depend on the new step or the case analysis above.

		if(function_current > function_best) {
			// Case a

			step_second = step_current;
			function_second = function_current;
			derivative_second = derivative_current;
		}
		else {
			if(opposite_sign_derivatives) {
				// Case c

				step_second = step_best;
				function_second = function_best;
				derivative_second = derivative_best;
			}

			// Case b/c

			step_best = step_current;
			function_best = function_current;
			derivative_best = derivative_current;
		}

		// Compute the new step and safeguard it

		step_trial = std::min(step_max, step_trial);
		step_trial = std::max(step_min, step_trial);
		step_current = step_trial;
		if(bracketed && bound) {
			// Mor\'e & Thuente (1994), p. 300, second formula

			// Note that the variable "delta" really has to be in
			// single precision to match the behavior of the original
			// Fortran implementation.
			const float delta = 0.66F /* < 1 */;
			const float c =
				step_best + delta * (step_second - step_best);

			step_current = (step_second > step_best) ?
				std::min(c, step_current) :
				std::max(c, step_current);
		}

		return STATUS_LINE_SEARCH_CONTINUE;
	}

	unsigned int reconstruction_filtering_t::solve_trust_region(
		float &newton_direction_pseudorapidity,
		float &newton_direction_azimuth,
		float &scaled_norm, bool &constrained,
		const float perp, const float gradient[], const float hessian[],
		const float determinant, const float gradient_norm,
		const float scaled_norm_limit, const bool constrain) const
	{
		I(scaled_norm_limit > 0);

		if(determinant > 0 || constrain) {
			// Definite Hessian
			solve_quadratic_unconstrained_2(
				newton_direction_pseudorapidity,
				newton_direction_azimuth, gradient, hessian,
				determinant);
		}
		else if(determinant < 0 && !constrain) {
			// Indefinite Hessian, handled by spectral decomposition
			if(gradient_norm < FLT_EPSILON * (1.0F + fabsf(perp))) {
				return STATUS_SADDLE_POINT;
				// FIXME: I need to handle saddle points, too
			}

			float modified_hessian[4] __attribute__ ((aligned(16)));

			modify_hessian_modified_cholesky(
				modified_hessian, hessian);

			const float modified_determinant =
				determinant_2(modified_hessian);

			// I(FEQ(1e-3 * modified_determinant, -1e-3 * determinant));

			solve_quadratic_unconstrained_2(
				newton_direction_pseudorapidity,
				newton_direction_azimuth, gradient, modified_hessian,
				modified_determinant);
		}
		else if(std::fpclassify(determinant) == FP_ZERO) {
			scaled_norm = NAN;
			return STATUS_HESSIAN_SINGULAR;
		}
		else if(!std::isfinite(determinant)) {
			scaled_norm = NAN;
			return STATUS_NOT_NUMERICAL;
		}
		else {
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": information: unhandled determinant = "
					  << determinant << std::endl;
		}

		const bool negative_semidefinite =
			determinant >= 0.0F && hessian[0] <= 0.0F;

		if(constrain || !negative_semidefinite) {
			static const float scale[2] = { 1.0F, 1.0F };

			scaled_norm = std::max(
				fabsf(scale[0] * newton_direction_pseudorapidity),
				fabsf(scale[1] * newton_direction_azimuth));

#if 0
			fprintf(stderr, "{{%s,%s},{%s,%s,%s}};\n",
					mathematica_form(gradient[0]),
					mathematica_form(gradient[1]),
					mathematica_form(hessian[0]),
					mathematica_form(hessian[1]),
					mathematica_form(hessian[2]));
			fprintf(stderr, "%s:%d: unconstrained: "
					"newton_direction = {%s,%s}; scaled_norm = %s; "
					"scaled_norm_limit = %s;\n",
					__FILE__, __LINE__,
					mathematica_form(newton_direction_pseudorapidity),
					mathematica_form(newton_direction_azimuth),
					mathematica_form(scaled_norm),
					mathematica_form(scaled_norm_limit));
#endif
			constrained = scaled_norm > scaled_norm_limit;
			if(constrained || !negative_semidefinite) {
				solve_quadratic_constrained_2(
					newton_direction_pseudorapidity,
					newton_direction_azimuth,
					-scaled_norm_limit / scale[0],
					scaled_norm_limit / scale[0],
					-scaled_norm_limit / scale[1],
					scaled_norm_limit / scale[1],
					gradient, hessian);
				scaled_norm = std::max(
					fabsf(scale[0] * newton_direction_pseudorapidity),
					fabsf(scale[1] * newton_direction_azimuth));
#if 0
				fprintf(stderr, "%s:%d: constrained: "
						"newton_direction = {%s,%s};\n",
						__FILE__, __LINE__,
						mathematica_form(newton_direction_pseudorapidity),
						mathematica_form(newton_direction_azimuth));
#endif
			}
		}
		else {
			scaled_norm = NAN;
		}

		return STATUS_CONTINUE;
	}

	unsigned int reconstruction_filtering_t::
	solve_line_search(
		float &step,
		const float newton_direction_pseudorapidity,
		const float newton_direction_azimuth,
		const jet_t &jet, const float perp, const float gradient[],
		const std::vector<track_t>::const_iterator &track_begin,
		const std::vector<track_t>::const_iterator &track_end,
		const collision_geometry_t &geometry) const
	{
		static const float function_tolerance = 1.0e-4F;
		static const float gradient_tolerance = 9.0e-1F;

		// Inexact line search

		float step_lower = 0.0F;	// best
		float step_upper = 0.0F;	// second
		float step_trial = 1.0F;	// current

		const float derivative =
			gradient[0] * newton_direction_pseudorapidity +
			gradient[1] * newton_direction_azimuth;
		const float function_test_derivative =
			-function_tolerance * derivative;
		const float derivative_test_derivative =
			gradient_tolerance * derivative;

		float perp_lower;
		float perp_upper;
		float perp_trial;
		float derivative_lower;
		float derivative_upper;
		float derivative_trial;

		// Initialization

		perp_lower = -perp;
		derivative_lower = -derivative;
		if(derivative_lower > 0) {
			std::cerr << __FILE__ << ':' << __LINE__ << ": error: "
				"ascending direction (derivative_lower = "
					  << mathematica_form(derivative_lower) << ";)"
					  << std::endl;
		}

		jet_t jet_upper = jet;
		float gradient_upper[2] __attribute__ ((aligned(8)));
		float hessian_upper[4] __attribute__ ((aligned(16)));

		jet_upper.momentum().pseudorapidity() +=
			step_upper * newton_direction_pseudorapidity;
		jet_upper.momentum().azimuth() = angular_range_reduce(
			jet_upper.momentum().azimuth() +
			step_upper * newton_direction_azimuth);
		evaluate_step(perp_upper, gradient_upper, hessian_upper,
					  jet_upper, track_begin, track_end,
					  geometry);
		perp_upper = -perp_upper;
		derivative_upper =
			-(gradient_upper[0] * newton_direction_pseudorapidity +
			  gradient_upper[1] * newton_direction_azimuth);
		if(derivative_upper > 0) {
			std::cerr << __FILE__ << ':' << __LINE__ << ": error: "
				"ascending direction (derivative_upper = "
					  << mathematica_form(derivative_upper) << ";)"
					  << std::endl;
		}

		static const int iteration_max = 64;
		float step_min = FLT_EPSILON;
		// The choice of step_max is chosen to be the typical maximum
		// meaningful step size for jet reconstruction, pi/4.
		float step_max = static_cast<float>(0.0625 * M_PI);
		float width[2] = {
			step_max - step_min, 2.0F * (step_max - step_min)
		};
		bool bracketed = false;
		unsigned int status = STATUS_LINE_SEARCH_CONTINUE;
		jet_t jet_trial;

		I(step_trial >= 0);

		for(int iteration = 0; iteration < iteration_max;
			iteration++) {
#if 0
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": position = {"
					  << mathematica_form(
							jet.momentum().pseudorapidity()) << ','
					  << mathematica_form(
							jet.momentum().azimuth()) << '}'
					  << std::endl;
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": newton_direction = {"
					  << mathematica_form(
							newton_direction_pseudorapidity) << ','
					  << mathematica_form(
							newton_direction_azimuth) << '}'
					  << std::endl;
#endif
			jet_trial = jet;
			jet_trial.momentum().pseudorapidity() +=
				step_trial * newton_direction_pseudorapidity;
			jet_trial.momentum().azimuth() = angular_range_reduce(
				jet_trial.momentum().azimuth() +
				step_trial * newton_direction_azimuth);
#if 0
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": position' = {"
					  << mathematica_form(
							jet_trial.momentum().pseudorapidity()) << ','
					  << mathematica_form(
							jet_trial.momentum().azimuth()) << '}'
					  << std::endl;
#endif

			float gradient_trial[2] __attribute__ ((aligned(8)));
			float hessian_trial[4] __attribute__ ((aligned(16)));

			evaluate_step(perp_trial, gradient_trial, hessian_trial,
						  jet_trial, track_begin, track_end,
						  geometry);
			perp_trial = -perp_trial;
			derivative_trial =
				-(gradient_trial[0] *
				  newton_direction_pseudorapidity +
				  gradient_trial[1] * newton_direction_azimuth);
#if 0
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": step = {"
					  << mathematica_form(step_lower) << ','
					  << mathematica_form(step_upper) << ','
					  << mathematica_form(step_trial) << '}'
					  << std::endl;
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": perp = {"
					  << mathematica_form(perp_lower) << ','
					  << mathematica_form(perp_upper) << ','
					  << mathematica_form(perp_trial) << '}'
					  << std::endl;
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": derivative = {"
					  << mathematica_form(derivative_lower) << ','
					  << mathematica_form(derivative_upper) << ','
					  << mathematica_form(derivative_trial) << '}'
					  << std::endl;
#endif

			// Test for convergence

			static const double x_tolerance = FLT_EPSILON;
			float function_test =
				-perp + step_trial * function_test_derivative;

			function_test *=
				(1.0F + (function_test >= 0.0F ? 1.0F : -1.0F) *
				 powf(FLT_EPSILON, 1.0F / 2.0F));
			if((perp_trial <= function_test &&
				fabsf(derivative_trial) <=
				derivative_test_derivative) ||
			   (bracketed &&
				(step_trial <= step_min || step_trial >= step_max ||
				 step_max - step_min < x_tolerance * step_max))) {
				status = STATUS_LINE_SEARCH_CONVERGENCE;
				break;
			}
			else if((step_trial >= step_max &&
					 perp_trial <= function_test &&
					 derivative <= function_test_derivative) ||
					(step_trial <= step_min &&
					 (perp_trial > function_test ||
					  derivative >= function_test_derivative))) {
#if 0
				std::cerr << __FILE__ << ':' << __LINE__
						  << ": step_trial = "
						  << mathematica_form(step_trial)
						  << "; step_min = "
						  << mathematica_form(step_min)
						  << "; step_max = "
						  << mathematica_form(step_max)
						  << "; perp_trial = "
						  << mathematica_form(perp_trial)
						  << "; function_test = "
						  << mathematica_form(function_test)
						  << "; derivative = "
						  << mathematica_form(derivative)
						  << "; function_test_derivative = "
						  << mathematica_form(
								function_test_derivative)
						  << ';' << std::endl;
#endif
				status = STATUS_LINE_SEARCH_STEP_LIMIT;
				break;
			}

			// Perform one line search iteration

#if 0
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": step = {"
					  << mathematica_form(step_lower) << ','
					  << mathematica_form(step_upper) << ','
					  << mathematica_form(step_trial) << '}'
					  << std::endl;
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": perp = {"
					  << mathematica_form(perp_lower) << ','
					  << mathematica_form(perp_upper) << ','
					  << mathematica_form(perp_trial) << '}'
					  << std::endl;
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": derivative = {"
					  << mathematica_form(derivative_lower) << ','
					  << mathematica_form(derivative_upper) << ','
					  << mathematica_form(derivative_trial) << '}'
					  << std::endl;
#endif

			const unsigned int step_status =
				line_search_more_thuente_step(
					step_lower, perp_lower, derivative_lower,
					step_upper, perp_upper, derivative_upper,
					step_trial, perp_trial, derivative_trial,
					bracketed, step_min, step_max);

			// If an unusual termination is to occur then let stp be
			// the lowest point obtained so far
			if(step_status != STATUS_LINE_SEARCH_CONTINUE) {
				step = step_lower;
				return step_status;
			}

			// Force a sufficient decrease in the size of the interval
			// of uncertainty
			if(bracketed) {
				const float step_upper_lower = step_upper - step_lower;
				const float fabs_step_upper_lower =
					fabsf(step_upper_lower);

				if(fabs_step_upper_lower >= 0.66F * width[1]) {
					step_trial = step_lower + 0.5F * step_upper_lower;
				}
				width[1] = width[0];
				width[0] = fabs_step_upper_lower;
				// Set the minimum and maximum steps to correspond to
				// the present interval of uncertainty
				if(step_lower < step_upper) {
					step_min = step_lower;
					step_max = step_upper;
				}
				else {
					step_min = step_upper;
					step_max = step_lower;
				}
			}
			else {
				if(step_trial < step_lower) {
					// Pathological case. There is no hope of
					// bracketing, take the constrained optimum and
					// continue.
					step = step_lower;
					return STATUS_LINE_SEARCH_STEP_LIMIT;
				}

				I(step_trial >= step_lower);

				static const float xtrapf = 4.0F;

				step_min = step_lower;
				step_max = step_trial +
					xtrapf * (step_trial - step_lower);
#if 1
				if(step_max < step_min) {
					std::cerr << __FILE__ << ':' << __LINE__
							  << ": information: step range adaption "
						"failed with the new step range = {"
							  << mathematica_form(step_min) << ','
							  << mathematica_form(step_max) << '}'
							  << std::endl;
				}
#endif

				I(step_max >= step_min);
			}
			// Force the step to be within bounds
			step_trial = std::min(
				std::max(step_trial, step_min), step_max);
#if 0
			fprintf(stderr, "%d: step = %f %f %f, %f %f, "
					"bracketed = %d\n",
					iteration, step_lower, step_upper,
					step_trial, step_min, step_max, bracketed);
#endif
		}
		step = step_trial;

		return status;
	}

	unsigned int reconstruction_filtering_t::refine_jet_maximum_step(
		jet_t &jet, snowmass_vector_t &previous_momentum,
		float &previous_determinant, float &scaled_norm_limit,
		const std::vector<track_t>::const_iterator &track_begin,
		const std::vector<track_t>::const_iterator &track_end,
		const float perp, const float gradient[],
		const float hessian[], const int iteration,
		const float precision, const collision_geometry_t &geometry,
		const bool biased_wolfe) const
	{
		// Termination condition
		const float gradient_norm = sqrtf(
			gradient[0] * gradient[0] + gradient[1] * gradient[1]);
		const float determinant =
			hessian[0] * hessian[1] - hessian[2] * hessian[2];

		if(iteration > 0) {
			// U1-U3 in Gill et al. (1981), section 8.2.3.2
			const float delta_perp = jet.momentum().perp() - perp;
			const float delta_position =
				jet.momentum().radial_distance(previous_momentum);
			const float accuracy_function =
				precision * (1.0F + fabsf(perp));
			const float accuracy_position = sqrtf(precision) *
				(1.0F + jet.momentum().radial_origin_distance());
			const float accuracy_gradient =
				cbrtf(precision) * (1.0F + fabsf(perp));

			if(delta_perp < accuracy_function &&
			   delta_position < accuracy_position &&
			   gradient_norm < accuracy_gradient) {
				return STATUS_CONVERGENCE;
			}

			// U4 in Gill et al. (1981), section 8.2.3.2, and with the
			// standard lower bound from section 8.5.1.3
			const float accuracy_gradient_norm =
				FLT_EPSILON * (1.0F + fabsf(perp));

			if(gradient_norm < accuracy_gradient_norm) {
				return STATUS_CONVERGENCE;
			}
#if 0
			// U5 in Gill et al. (1981), section 8.2.3.2
			if(previous_determinant > 0 && determinant < 0) {
				return STATUS_HESSIAN_INDEFINITE;
			}
#endif
		}
		previous_momentum = jet.momentum();
		previous_determinant = determinant;

#if 0
		std::cerr << __FILE__ << ':' << __LINE__
				  << ": information: iteration = " << iteration
				  << "; jet.momentum = {"
				  << mathematica_form(jet.momentum()[1]) << ','
				  << mathematica_form(jet.momentum()[2]) << ','
				  << mathematica_form(jet.momentum()[3] / M_PI)
				  << " Pi};" << std::endl;
#endif

		float newton_direction_pseudorapidity;
		float newton_direction_azimuth;
		float scaled_norm;
		bool constrained;

		I(scaled_norm_limit > 0);

		// Jet reconstruction has a natural scale, beyond which
		// constrained steps does not make sense
		scaled_norm_limit =
			std::max(static_cast<float>(0.0625 * M_PI),
					 scaled_norm_limit);

#if 0
		const int trust_region_status =
#endif
			solve_trust_region(
				newton_direction_pseudorapidity,
				newton_direction_azimuth, scaled_norm, constrained,
				perp, gradient, hessian, determinant, gradient_norm,
				scaled_norm_limit, biased_wolfe);
		if(biased_wolfe) {
			const float q = evaluate_quadratic(
				newton_direction_pseudorapidity,
				newton_direction_azimuth, gradient, hessian);

			jet_t jet_step = jet;

			jet_step.momentum().pseudorapidity() +=
				newton_direction_pseudorapidity;
			jet_step.momentum().azimuth() = angular_range_reduce(
				jet_step.momentum().azimuth() +
				newton_direction_azimuth);

			float perp_step;
			float gradient_step[2] __attribute__ ((aligned(8)));
			float hessian_step[4] __attribute__ ((aligned(16)));

			evaluate_step(perp_step, gradient_step, hessian_step,
						  jet_step, track_begin, track_end,
						  geometry);

			const float r = (perp_step - perp) / q;
			float step;

#if 0
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": iteration = " << iteration
					  << "; q = " << mathematica_form(q)
					  << "; r = " << mathematica_form(r)
					  << ';' << std::endl;
#endif

			const int line_search_status = solve_line_search(
				step, newton_direction_pseudorapidity,
				newton_direction_azimuth, jet, perp, gradient,
				track_begin, track_end, geometry);

			if(line_search_status !=
			   STATUS_LINE_SEARCH_CONVERGENCE &&
			   line_search_status !=
			   STATUS_LINE_SEARCH_STEP_LIMIT) {
#if 0
				std::cerr << __FILE__ << ':' << __LINE__
						  << ": information: status = "
						  << status_str(line_search_status)
						  << std::endl;
#endif
				if(line_search_status ==
				   STATUS_LINE_SEARCH_NOT_NUMERICAL) {
					return STATUS_NOT_NUMERICAL;
				}
			}

			static const float trust_region_threshold_upper = 0.9;
			static const float nu = 1.0F;
			static const float trust_region_increase = 2.0F;

			I(scaled_norm_limit > 0);

			if(r > trust_region_threshold_upper) {
				// ||N_k s_k|| >= (1 / nu) delta_k unless s_k =
				// -B_k^(-1) g
				const float nu_hat = 0.5F * (1.0F + nu);

				scaled_norm_limit = std::max(
					scaled_norm_limit,
					std::max(step * nu_hat * scaled_norm,
							 trust_region_increase * scaled_norm));

				I(scaled_norm_limit > 0);
			}
			else {
				IG(scaled_norm >= 0,
				   scaled_norm * nu <= scaled_norm_limit);

				// As a modification to the original biased Wolfe
				// algorith, skip the update if the line search fails
				// with step <= 0
				if(step > 0) {
					scaled_norm_limit = step *
						(scaled_norm * nu <= scaled_norm_limit ?
						 scaled_norm : scaled_norm_limit);
				}

				IG(scaled_norm_limit > 0, scaled_norm > 0);
			}
		}

#if 0
		std::cerr << __FILE__ << ':' << __LINE__
				  << ": information: iteration = " << iteration
				  << "; gradient = {"
				  << mathematica_form(gradient[0]) << ','
				  << mathematica_form(gradient[1])
				  << "}; hessian = {" << mathematica_form(hessian[0])
				  << ',' << mathematica_form(hessian[1]) << ','
				  << mathematica_form(hessian[2]) << "};"
				  << std::endl;
		std::cerr << __FILE__ << ':' << __LINE__
				  << ": information: iteration = " << iteration
				  << "; jet.momentum = {"
				  << mathematica_form(jet.momentum()[2]) << ','
				  << mathematica_form(jet.momentum()[3] / M_PI)
				  << " Pi}; newton_direction = {"
				  << mathematica_form(
						newton_direction_pseudorapidity) << ','
				  << mathematica_form(
						newton_direction_azimuth / M_PI) << " Pi};"
				  << std::endl;
#endif
#if 0
		// Acceptance limit of the background model. This has to been
		// applied before trust-region and line search methods.
		if(geometry.in_range()) {
			constrain_newton_direction_pseudorapidity(
				newton_direction_pseudorapidity,
				newton_direction_azimuth, jet.momentum()[2]);
		}
#endif

		// Apply the Newton step
		jet.momentum()[1] = perp;
		jet.momentum()[2] += newton_direction_pseudorapidity;
		jet.momentum()[3] = angular_range_reduce(
			jet.momentum()[3] + newton_direction_azimuth);

		return STATUS_CONTINUE;
	}

	void reconstruction_filtering_t::
	refine_jet_maximum(
		jet_t &jet,
		const std::vector<track_t>::const_iterator &track_begin,
		const std::vector<track_t>::const_iterator &track_end,
		const collision_geometry_t &geometry) const
	{
		static const float precision = 4.0F * FLT_EPSILON;
		static const int max_iteration = 64;
		snowmass_vector_t previous_momentum;
		float previous_determinant;
		float perp;
		float gradient[2] __attribute__ ((aligned(8)));
		float hessian[4] __attribute__ ((aligned(16)));
		int iteration;
		int status = STATUS_CONVERGENCE;
		const jet_t unrefined_jet = jet;

		float scaled_norm_limit =
			2.0F * std::pow(_standard_deviation, 4);

		for(iteration = 0; iteration < max_iteration; iteration++) {
			evaluate_step(perp, gradient, hessian, jet,
						  track_begin, track_end, geometry);
#if 0
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": scale = "
					  << mathematica_form(
							_background_perp->
							vertex_centrality_dependence(geometry))
					  << "; x = {"
					  << mathematica_form(
							jet.momentum().pseudorapidity())
					  << ','
					  << mathematica_form(
							jet.momentum().azimuth())
					  << "}; gradient = {"
					  << mathematica_form(gradient[0]) << ','
					  << mathematica_form(gradient[1])
					  << "}; hessian = {"
					  << mathematica_form(hessian[0]) << ','
					  << mathematica_form(hessian[1]) << ','
					  << mathematica_form(hessian[2]) << '}'
					  << std::endl;
#endif
			status = refine_jet_maximum_step(
				jet, previous_momentum, previous_determinant,
				scaled_norm_limit, track_begin, track_end,
				perp, gradient, hessian, iteration, precision,
				geometry, true);
			if(status != STATUS_CONTINUE)
				break;
		}

		if(status == STATUS_CONVERGENCE) {
			jet.momentum().set_lightlike_time();
#if 0
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": information: distance = {"
					  << mathematica_form(
							jet.momentum().pseudorapidity() -
							jet.seed_momentum().pseudorapidity())
					  << ','
					  << mathematica_form(angular_range_reduce(
							jet.momentum().azimuth() -
							jet.seed_momentum().azimuth()))
					  << '}' << std::endl;
#endif
		}
		else {
			// Activate the following to test for mismatch between
			// discrete filtering and continuous maximization, Cu + Cu
			// events for instance will naturally produce many
			// signular Hessians.
#if 0
			std::cerr << __FILE__ << ':' << __LINE__
					  << ": information: status = " << status_str(status)
					  << "; position = {"
					  << mathematica_form(
							unrefined_jet.momentum().perp())
					  << ','
					  << mathematica_form(
							unrefined_jet.momentum().
							pseudorapidity())
					  << ','
					  << mathematica_form(
							unrefined_jet.momentum().azimuth())
					  << "};"
					  << std::endl;
#endif
			jet.momentum().perp() = 0;
		}
	}

}
