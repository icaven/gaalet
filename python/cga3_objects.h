//
// Conversions between CGA3 rotors, translators and motors to/from rotation and translation vectors
//

#ifndef CGA3_OBJECTS_H
#define CGA3_OBJECTS_H

#include <vector>
#include "cga3d.h"


template <class ELEMENT>
typename cga3<ELEMENT>::Rotor rotor_from_rotation_vector(const std::vector<ELEMENT> &rvec)
{
	ELEMENT r1 = rvec[0];
	ELEMENT r2 = rvec[1];
	ELEMENT r3 = rvec[2];

	const auto& layout = cga3<ELEMENT>::get();

	ELEMENT angle_squared = r1 * r1 + r2 * r2 + r3 * r3;
	ELEMENT theta = sqrt(angle_squared);

	// Arbitrarily initialize the axis to be the z-axis in case the angle is 0
	typename cga3<ELEMENT>::Rotor axis = layout.e12;
	if (theta > 0.0)
	{
		axis = r3 * layout.e12 + r2 * layout.e31 + r1 * layout.e23;
		axis = 1./ theta * axis;
	}
	else
	{
		theta = ELEMENT(0.0);   // Ensure that the derivatives are zero
	}

	typename cga3<ELEMENT>::Scalar c = cos(0.5 * theta);
	typename cga3<ELEMENT>::Rotor s = sin(0.5 * theta) * axis;
	typename cga3<ELEMENT>::Rotor rotor = c - s;

	return rotor;
}


template <class ELEMENT, class E>
std::vector<ELEMENT> rotation_vector_from_rotor(const gaalet::expression<E>& r)
{
	const auto& layout = cga3<ELEMENT>::get();
	std::vector<ELEMENT> rvec = {ELEMENT(0), ELEMENT(0), ELEMENT(0)};
	if (::magnitude2(::grade<2>(r)).template element<0>() > 0 || ::grade<0>(r).template element<0>() != 1)
	{
		ELEMENT half_angle = acos(::grade<0>(r).template element<0>());
		typename cga3<ELEMENT>::Rotor axis = eval(-1 * layout.e12);    // Set an arbitrary axis if the rotation angle is 0
		if (half_angle != ELEMENT(0))
		{
			axis = eval(1. / sin(half_angle) * ::grade<2>(r));
		}
		typename cga3<ELEMENT>::EuclideanVector euc_vector = ::grade<1>(axis * layout.I3);
		rvec[0] = 2. * half_angle * euc_vector.template element<vector_conf(1)>();
		rvec[1] = 2. * half_angle * euc_vector.template element<vector_conf(2)>();
		rvec[2] = 2. * half_angle * euc_vector.template element<vector_conf(3)>();
	}
	return rvec;
}

template <class ELEMENT>
typename cga3<ELEMENT>::Translator translator_from_translation_vector(const std::vector<ELEMENT>& tvec)
{
	const auto& layout = cga3<ELEMENT>::get();
	typename cga3<ELEMENT>::Translator translator;
	translator = layout.one - 0.5 * (tvec[0] * layout.e1inf + tvec[1] * layout.e2inf + tvec[2] * layout.e3inf);
	return translator;
}

template <class ELEMENT, class E>
typename cga3<ELEMENT>::Translator generate_translation_rotor(const gaalet::expression<E>& t)
{
	const auto& layout = cga3<ELEMENT>::get();
	typename cga3<ELEMENT>::Translator translator;

	translator = layout.one + 0.5 * layout.einf * t;
	return translator;
}


template <class ELEMENT>
typename cga3<ELEMENT>::Motor motor_from_rt_parameters(const std::vector<ELEMENT> &extrinsic_parameters)
{
	ELEMENT t1 = extrinsic_parameters[3];
	ELEMENT t2 = extrinsic_parameters[4];
	ELEMENT t3 = extrinsic_parameters[5];

	std::vector<ELEMENT> translation_components({t1, t2, t3});
	typename cga3<ELEMENT>::Translator translator = translator_from_translation_vector<ELEMENT>(translation_components);
	typename cga3<ELEMENT>::Rotor rotor = rotor_from_rotation_vector(extrinsic_parameters);

	return translator * rotor;
}

template <class ELEMENT, class E>
std::vector<ELEMENT> rt_parameters_from_motor(const gaalet::expression<E>& motor)
{
	const auto& layout = cga3<ELEMENT>::get();
	std::vector<ELEMENT> rt_parameters {ELEMENT(0), ELEMENT(0), ELEMENT(0), ELEMENT(0), ELEMENT(0), ELEMENT(0)};

	// First, decompose the motor into a pure rotation rotor and a component representing the translation
	// p. 383-384 of Dorst, Fontijne, and Mann
	typename cga3<ELEMENT>::Rotor r = -1 * layout.eo & (motor * layout.einf);   // The pure rotation rotor

	if (magnitude2(::grade<2>(r)).template element<0>() == ELEMENT(0) && ::grade<0>(r).template element<0>() == ELEMENT(-1))
	{
		// No unique logarithm - return the parameters for the identity transform
		return rt_parameters;
	}

	auto t = -2 * (layout.eo & motor) * inverse(r);   // This can be the sum of a vector and a trivector
	if (magnitude2(::grade<2>(r)).template element<0>() == ELEMENT(0) && ::grade<0>(r).template element<0>() == ELEMENT(1))
	{
		// The rotation vector is zero, so just convert the translator
		rt_parameters[3] = ::grade<1>(t).template element<vector_conf(1)>();
		rt_parameters[4] = ::grade<1>(t).template element<vector_conf(2)>();
		rt_parameters[5] = ::grade<1>(t).template element<vector_conf(3)>();
		return rt_parameters;
	}

	auto rotation_vector = rotation_vector_from_rotor<ELEMENT>(r);
	rt_parameters[0] = rotation_vector[0];
	rt_parameters[1] = rotation_vector[1];
	rt_parameters[2] = rotation_vector[2];

	// Even though this decomposition works, we still don't have a vector describing the translation
	// since t can be the sum of a vector and a trivector
	// Determine the translation vector tvec and the rotor for the screw motion
	// such that T * r == motor, where T is the translator from the tvec, and r is the pure rotor
	auto plane_of_rotation = cga3<ELEMENT>::normalise(r);        // The plane of rotation
	auto w = (t ^ plane_of_rotation) * inverse(plane_of_rotation);
	auto u = (t & plane_of_rotation) * inverse(plane_of_rotation);
	auto v = (1/(1 - (r*r).template element<0>())) * u;
	typename cga3<ELEMENT>::Translator t_v = generate_translation_rotor<ELEMENT>(v);
	typename cga3<ELEMENT>::Translator t_r = generate_translation_rotor<ELEMENT>(w) * t_v * r * ~t_v;
	typename cga3<ELEMENT>::Vector translation_vector = -2. * (layout.eo & (t_r * ~r));

	rt_parameters[3] = ::grade<1>(translation_vector).template element<vector_conf(1)>();
	rt_parameters[4] = ::grade<1>(translation_vector).template element<vector_conf(2)>();
	rt_parameters[5] = ::grade<1>(translation_vector).template element<vector_conf(3)>();
	return rt_parameters;
}


template <class ELEMENT, class L1, class L2>
ELEMENT angle_between_line_pair(const gaalet::expression<L1>& line1, const gaalet::expression<L2>& line2)
{
	typename cga3<ELEMENT>::Rotor rotor_squared = line1 * inverse(line2);
	ELEMENT scalar_part = ::grade<0>(rotor_squared).template element<0>();
	ELEMENT magnitude_bivector_part = magnitude(::grade<2>(rotor_squared)).template element<0>();

	// This will return an angles in the range 0 to pi
	return atan2(magnitude_bivector_part, scalar_part);
}

#endif //CGA3_OBJECTS_H
