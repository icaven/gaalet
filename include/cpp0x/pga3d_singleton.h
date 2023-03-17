//
// Class to encapsulate operators and objects of 3d Plane based Geometric Algebra
// This allows the algebra to be used with different element types at the same time.
//

#ifndef PGA3D_H
#define PGA3D_H

#pragma once

#include <cstdarg>
#include <cassert>
#include <vector>
#include <opencv2/core/matx.hpp>

#include "gaalet.h"
#include "pga3.h"
#include "pga3_dual.h"
#include "pga3_norm.h"
#include "pga3_exponential.h"
#include "pga3_logarithm.h"
#include "magnitude.h"

/// @brief Determine if two elements are close in value
/// @param a
/// @param b
/// @param rtol
/// @param atol
/// @return
template <class ELEMENT>
inline bool isclose(ELEMENT a, ELEMENT b, ELEMENT rtol = ELEMENT(1e-5), ELEMENT atol = ELEMENT(DBL_EPSILON))
{
	return abs(a - b) <= (atol + rtol * abs(b));
}


/// Define some configuration lists for common multivector types
#define DIRECTION_CONF \
	blade_conf(0, 3, 2), blade_conf(0, 1, 3), blade_conf(0, 2, 1)
#define POINT_CONF \
	DIRECTION_CONF, blade_conf(1, 2, 3)
#define EUC_VECTOR_CONF \
	pga3::vector_conf(1), pga3::vector_conf(2), pga3::vector_conf(3)
#define VECTOR_CONF \
  pga3::vector_conf(1), pga3::vector_conf(2), pga3::vector_conf(3), pga3::vector_conf(0)
#define BIVECTOR_CONF \
	blade_conf(0, 1),  blade_conf(0, 2), blade_conf(0, 3), \
	blade_conf(2, 3), blade_conf(1, 3), blade_conf(1, 2)
#define ROTOR_CONF \
	0, blade_conf(2, 3), blade_conf(1, 3), blade_conf(1, 2)
#define TRANSLATOR_CONF \
	0, blade_conf(0, 1),  blade_conf(0, 2), blade_conf(0, 3)
#define MOTOR_CONF \
	0, BIVECTOR_CONF, pga3::pseudoscalar_conf
#define PLANE_CONF VECTOR_CONF
#define LINE_CONF BIVECTOR_CONF


template <typename ELEMENT>
class pga3d {

public:

	///@brief Definition of the 4 dimensional plane based geometric algebra space
	typedef gaalet::algebra<gaalet::signature<3, 0, 1>, ELEMENT> algebra;
	static constexpr gaalet::conf_t pseudoscalar_conf = pga3::pseudoscalar_conf;

	// The multivector type for the algebra
	template<gaalet::conf_t... elements>
	struct mv {
		typedef typename algebra::template mv<elements...>::type type;
	};

	/// Frequently used multivector types
	typedef typename mv<0x0>::type Scalar;
	typedef typename mv<POINT_CONF>::type Point;
	typedef typename mv<EUC_VECTOR_CONF>::type EuclideanVector;
	typedef typename mv<DIRECTION_CONF>::type Direction;
	typedef typename mv<BIVECTOR_CONF>::type Bivector;
	typedef typename mv<LINE_CONF>::type Line;
	typedef typename mv<PLANE_CONF>::type Plane;
	typedef typename mv<ROTOR_CONF>::type Rotor;
	typedef typename mv<TRANSLATOR_CONF>::type Translator;
	typedef typename mv<MOTOR_CONF>::type Motor;

	/// The basis elements of the algebra

	/// Scalar
	Scalar one;

	/// Vectors
	typename mv<pga3::vector_conf(0)>::type e0;
	typename mv<pga3::vector_conf(1)>::type e1;
	typename mv<pga3::vector_conf(2)>::type e2;
	typename mv<pga3::vector_conf(3)>::type e3;

	/// Bivectors used when representing motors
	typename mv<blade_conf(1, 2)>::type e12;    // rotation around z-axis
	typename mv<blade_conf(1, 3)>::type e31;    // rotation around y-axis
	typename mv<blade_conf(2, 3)>::type e23;    // rotation around x-axis
	typename mv<blade_conf(0, 1)>::type e01;    // translation in x-axis direction
	typename mv<blade_conf(0, 2)>::type e02;    // translation in y-axis direction
	typename mv<blade_conf(0, 3)>::type e03;    // translation in z-axis direction

	/// Other names for the bivectors
	typename mv<blade_conf(1, 2)>::type i = e12;      // z-axis      - see Section 7.4.3 of [Gunn2011]
	typename mv<blade_conf(1, 3)>::type j = e31;      // y-axis
	typename mv<blade_conf(2, 3)>::type k = e23;      // x-axis
	typename mv<blade_conf(0, 1)>::type dual_k = e01;
	typename mv<blade_conf(0, 2)>::type dual_j = e02;
	typename mv<blade_conf(0, 3)>::type dual_i = e03;

	/// Name of trivectors
	typename mv<blade_conf(1, 2, 3)>::type E0;
	typename mv<blade_conf(0, 3, 2)>::type E1;
	typename mv<blade_conf(0, 1, 3)>::type E2;
	typename mv<blade_conf(0, 2, 1)>::type E3;

	/// Pseudoscalar for Euclidean space
	typename mv<blade_conf(1, 2, 3)>::type I3;

	/// The pseudoscalar (representing all space)
	typename mv<pseudoscalar_conf>::type I;

	/// Even grade configuration list
	// The output configuration list is composed of all the bivector elements and the scalar and pseudoscalar
	typedef typename insert_element<0,
					typename insert_element<pga3::I_CONF,
					typename insert_element<pga3::J_CONF,
					typename insert_element<pga3::K_CONF,
					typename insert_element<pga3::DUAL_K_CONF,
					typename insert_element<pga3::DUAL_J_CONF,
					typename insert_element<pga3::DUAL_I_CONF,
					typename insert_element<pseudoscalar_conf,
					cl_null>::clist>::clist>::clist>::clist>::clist>::clist>::clist>::clist even_grade_clist;

	/// Vector configuration list
	typedef typename insert_element<pga3::vector_conf(0),
					typename insert_element<pga3::vector_conf(1),
					typename insert_element<pga3::vector_conf(2),
					typename insert_element<pga3::vector_conf(3), 
					cl_null>::clist>::clist>::clist>::clist vector_clist;

	/// Bivector configuration list
	typedef typename insert_element<pga3::I_CONF,
					typename insert_element<pga3::J_CONF,
					typename insert_element<pga3::K_CONF,
					typename insert_element<pga3::DUAL_K_CONF,
					typename insert_element<pga3::DUAL_J_CONF,
					typename insert_element<pga3::DUAL_I_CONF,
					cl_null>::clist>::clist>::clist>::clist>::clist>::clist bivector_clist;


	/// Constructor - initialises the basis elements
	pga3d()
	{
		one = {ELEMENT(1.0)};
		e0 = {ELEMENT(1.0)};
		e1 = {ELEMENT(1.0)};
		e2 = {ELEMENT(1.0)};
		e3 = {ELEMENT(1.0)};
		i = e12 = {eval(e1 ^ e2)};
		j = e31 = {eval(e3 ^ e1)};
		k = e23 = {eval(e2 ^ e3)};
		dual_k = e01 = {eval(e0 ^ e1)};
		dual_j = e02 = {eval(e0 ^ e2)};
		dual_i = e03 = {eval(e0 ^ e3)};
		E0 = {eval(e1 ^ e2 ^ e3)};      // The origin point
		E1 = {eval(e0 ^ e3 ^ e2)};
		E2 = {eval(e0 ^ e1 ^ e3)};
		E3 = {eval(e0 ^ e2 ^ e1)};
		I3 = {eval(e1 ^ e2 ^ e3)};
		I = {eval(e0 ^ e1 ^ e2 ^ e3)};
	}

	/// No copy, move, or assignment operator since this is a singleton
	pga3d(const pga3d& other) = delete;
	pga3d(const pga3d&& other) = delete;
	pga3d& operator=(const pga3d& other) = delete;

	/// Utility functions
	template <class L, class R>
	bool isclose(const gaalet::expression<L>& l, const gaalet::expression<R>& r,
							 ELEMENT rtol = ELEMENT(1e-5), ELEMENT atol = ELEMENT(DBL_EPSILON)) const
	{
		auto n = norm(l - r);
		return ::isclose(n, ELEMENT(0), rtol, atol);
	}

	inline auto make_motor(const ELEMENT &scalar, const ELEMENT &dual_k_value, const ELEMENT &dual_j_value,
	                       const ELEMENT &dual_i_value, const ELEMENT &i_value, const ELEMENT &j_value,
	                       const ELEMENT &k_value, const ELEMENT &pseudoscalar_value) const
	{
		return scalar * one
		       + dual_k_value * dual_k
		       + dual_j_value * dual_j
		       + dual_i_value * dual_i
		       + i_value * i
		       + j_value * j
		       + k_value * k
		       + pseudoscalar_value * I;
	}

	/// Return the Study number and the axis (a normalized bivector) associated with the given bivector
	// See: Section 7.7 of Gunn, Charles, "Geometry, Kinematics, and Rigid Body Mechanics in Cayley-Klein Geometries" [Gunn2011]
	// Normalize the bivector by multiplying it by the inverse of the square root of the associated Study number z, and the result
	// is B, which is the axis of the given bivector
	// The paper by de Keninck and Roelfs provides a similar algorithm
	std::pair<std::pair<ELEMENT, ELEMENT>, Bivector>
	bivector_axis(const Bivector& bivector) const {
		auto z = ::eval(bivector * (~bivector));
		ELEMENT a = ::grade<0>(z).template element<0>();
		auto axis = ::eval(
						bivector);   // In the case that the inverse sqrt of z doesn't exist (when a == 0), the B is equal to the bivector
		ELEMENT b = ELEMENT(0);
		std::pair<ELEMENT, ELEMENT> study_number = {ELEMENT(0), ELEMENT(0)};
		if (!::isclose(a, ELEMENT(0))) {
			auto sqrt_a = sqrt(a);
			b = ::grade<4>(z).template element<pseudoscalar_conf>();
			auto inv_sqrt_z = sqrt_a / a * one - ((b) / (ELEMENT(2) * a * sqrt_a)) * I;
			axis = inv_sqrt_z * axis;
			study_number = std::make_pair(sqrt_a, b / (ELEMENT(2) * sqrt_a));
			return std::make_pair(study_number, axis);
		}
		return std::make_pair(study_number, axis);
	}

	Point make_point(const ELEMENT& x, const ELEMENT& y, const ELEMENT& z) const
	{
		return E0 + x * E1 + y * E2 + z * E3;
	}


	/// @brief Lines can be defined by Plücker coordinates.
	inline
	auto make_line(ELEMENT px, ELEMENT py, ELEMENT pz, ELEMENT dx, ELEMENT dy, ELEMENT dz) const
	{
		return dx * dual_k + dy * dual_j + dz * dual_i + px * k + py * j + pz * i;
	}

	///@brief Raise a Euclidean point to a pga3 point
	Point up(const std::vector<ELEMENT>& euc_pt) const
	{
		return E0 + ELEMENT(euc_pt[0]) * E1 + ELEMENT(euc_pt[1]) * E2 + ELEMENT(euc_pt[2]) * E3;
	}

	Direction direction(const std::vector<ELEMENT>& euc_pt) const
	{
		return ELEMENT(euc_pt[0]) * E1 + ELEMENT(euc_pt[1]) * E2 + ELEMENT(euc_pt[2]) * E3;
	}

	Point up(const cv::Vec3d& euc_pt) const
	{
		return E0 + ELEMENT(euc_pt[0]) * E1 + ELEMENT(euc_pt[1]) * E2 + ELEMENT(euc_pt[2]) * E3;
	}

	Point up(const cv::Vec4d& euc_pt) const
	{
		return E0 + ELEMENT(euc_pt[0]) * E1 + ELEMENT(euc_pt[1]) * E2 + ELEMENT(euc_pt[2]) * E3;
	}

	template <class E>
	auto dual(const E& x) const
	{
		return pga3::dual(x);
	}

	// Implementation for exp(), log(), sqrt(), based on reference implementations in:
	// de Keninck, Steven and Roelfs, Martin "Normalization, Square Roots, and the Exponential and Logarithmic
	//                                        Maps in Geometric Algebras of Less than 6D"


	template <class E>
	pga3d<ELEMENT>::Motor bivector_exp(const gaalet::expression<E>& x) const
	{
		auto eval_x = eval(x);
		// Extract the components of the rotor
		auto i_value = eval_x.template element<pga3::I_CONF>();
		// Need to negate the j value since the value is stored for even permutation of the basis vectors,
		// but it is for an odd permutation.  This is an implementation detail.
		auto j_value = ELEMENT(-1) * eval_x.template element<pga3::J_CONF>();
		auto k_value = eval_x.template element<pga3::K_CONF>();
		auto dual_i_value = eval_x.template element<pga3::DUAL_I_CONF>();
		auto dual_j_value = eval_x.template element<pga3::DUAL_J_CONF>();
		auto dual_k_value = eval_x.template element<pga3::DUAL_K_CONF>();

		auto l = i_value * i_value + j_value * j_value + k_value * k_value;
		if (l == ELEMENT(0))
		{
			return one + dual_k_value * dual_k + dual_j_value * dual_j + dual_i_value * dual_i;
		}
		auto m = dual_k_value * k_value + dual_j_value * j_value + dual_i_value * i_value;
		auto a = sqrt(l);
		auto c = cos(a);
		auto s = sin(a) / a;
		auto t = m / l * (c - s);

		auto exp_x = make_motor(c,
														s * dual_k_value + t * k_value, 
														s * dual_j_value + t * j_value, 
														s * dual_i_value + t * i_value,
														s * i_value, 
														s * j_value, 
														s * k_value,
														m * s);
		return exp_x;

	}

	/// Objects to select the correct exponential function at compile time using the grade

	template<class A, int ET = pga3::detail::exponential_evaluation_type<A>::value>
	struct exponential : public expression<exponential<A>>
	{
		static_assert(ET!=-1, "no method for evaluating this type of multivector implemented");
	};

	// PGA3 exponential:bivectors
	template <class A>
	struct exponential<A, 2> : public expression<exponential<A>>
	{
		static const conf_t pseudoscalar_conf = Power<2, A::metric::dimension>::value - 1;

		// The output configuration list is composed of all of the bivector elements and the scalar and pseudoscalar
		typedef even_grade_clist clist;
		typedef typename A::metric metric;
		typedef typename A::element_t element_t;

		explicit exponential(const A &a_)
						: a(a_)
		{
		}

		auto operator()()
		{
			auto& layout = pga3d<ELEMENT>::get();
			return layout.bivector_exp(a);
		}
	private:
		A a;
	};

	// Scalar exponential
	template <class A>
	struct exponential<A, 0> : public expression<exponential<A>>
	{
		// The output configuration list is composed of the scalar grade only
		typedef typename insert_element<0, cl_null>::clist clist;

		typedef typename A::metric metric;

		typedef typename A::element_t element_t;

		explicit exponential(const A &a_)
						: a(a_)
		{
		}

		auto operator()()
		{
			return ::exp(a);
		}
	private:
		A a;
	};

	template <class E>
	auto exp(const gaalet::expression<E>& x) const
	{
		return exponential<E>(x)();
	}

	/// @brief Logarithm of an even grade multivector
	/// @tparam E 	The class of the value
	/// @param x 		The value
	/// @return 	The logarithm bivector
	template <class E>
	pga3d<ELEMENT>::Bivector log_of_rotor(const gaalet::expression<E>& x) const
	{
		auto bivector_part = ::grade<2>(x);
		auto scalar_part = ::grade<0>(x).template element<0>();
		auto dual_i_value = bivector_part.template element<pga3::DUAL_I_CONF>();
		auto dual_j_value = bivector_part.template element<pga3::DUAL_J_CONF>();
		auto dual_k_value = bivector_part.template element<pga3::DUAL_K_CONF>();
		if (scalar_part == ELEMENT(1))
		{
			// The multivector is a translator
			return dual_k_value * dual_k + dual_j_value * dual_j + dual_i_value * dual_i;
		}
		// The multivector is a motor
		auto pseudoscalar_part = ::grade<4>(x).template element<pseudoscalar_conf>();

		ELEMENT a = ELEMENT(1) / (ELEMENT(1) - scalar_part * scalar_part);
		ELEMENT b = acos(scalar_part) * sqrt(a);
		ELEMENT c = a * pseudoscalar_part * (ELEMENT(1) - scalar_part * b);

		// Extract the components of the rotor
		auto i_value = bivector_part.template element<pga3::I_CONF>();
		// Need to negate the j value since the value is stored for even permutation of the basis vectors,
		// but it is for an odd permutation.  This is an implementation detail.
		auto j_value = ELEMENT(-1) * bivector_part.template element<pga3::J_CONF>();
		auto k_value = bivector_part.template element<pga3::K_CONF>();

		auto log_x = (c * k_value + b * dual_k_value) * dual_k + (c * j_value + b * dual_j_value) * dual_j +
		             (c * i_value + b * dual_i_value) * dual_i + b * i_value * i + b * j_value * j + b * k_value * k;
		return log_x;
	}

	/// Objects to select the correct logarithm function at compile time using the grade
	template<class A, int ET = pga3::detail::logarithm_evaluation_type<A>::value>
	struct logarithm : public expression<logarithm<A>>
	{
		static_assert(ET!=-1, "no method for evaluating this type of multivector implemented");
	};

	/// PGA3 logarithm of even grade multivectors
	template <class A>
	struct logarithm<A, 1> : public expression<logarithm<A>>
	{
		// The output configuration list is composed of all of the bivector elements
		typedef bivector_clist clist;
		typedef typename A::metric metric;
		typedef typename A::element_t element_t;

		explicit logarithm(const A &a_)
						: a(a_)
		{
		}

		auto operator()()
		{
			auto& layout = pga3d<ELEMENT>::get();
			return layout.log_of_rotor(a);
		}
	private:
		A a;
	};

	/// Scalar logarithm
	template <class A>
	struct logarithm<A, 0> : public expression<logarithm<A>>
	{
		// The output configuration list is composed of the scalar grade only
		typedef typename insert_element<0, cl_null>::clist clist;

		typedef typename A::metric metric;

		typedef typename A::element_t element_t;

		explicit logarithm(const A &a_)
						: a(a_)
		{
		}

		auto operator()()
		{
			return ::log(a);
		}
	private:
		A a;
	};

	/// @brief The logarithm of a PGA3 multivector (only even grade multivectors and scalars supported)
	/// @tparam E 	The expression class of the value
	/// @param x 		The expression
	/// @return 	The logarithm
	template <class E>
	auto log(const gaalet::expression<E>& x) const
	{
		return logarithm<E>(x)();
	}

	/// Objects to select the correct normalize function at compile time using the grade
	// check for multivectors of that can be normalized
	// value=1 : vector grade only
	// value=2 : even grade only
	// value=0 : scalar, vector, trivector, pseudoscalar or a mix of grades (sometimes this is chosen because of 
	//           coefficients that are close to 0)
	template<class A>
	struct normalize_evaluation_type {
		static const int value =
						check_vector<typename A::clist>::value ? 1
						: check_even_grade<typename A::clist>::value ? 2
						: 0;
	};

	// This object should never be created, but it may be used if the normalize_evaluation_type is changed to handle
	// multivectors differently
	template<class A, int ET = normalize_evaluation_type<A>::value>
	struct normalize_mv : public expression<normalize_mv<A>>
	{
		static_assert(ET!=-1, "no method for evaluating this type of multivector implemented");
	};

	// Hyperplane normalization
	template <class A>
	struct normalize_mv<A, 1> : public expression<normalize_mv<A>>
	{
		// The output configuration list is composed of all of the vector elements
		typedef typename insert_element<pga3::vector_conf(0),
						typename insert_element<pga3::vector_conf(1),
										typename insert_element<pga3::vector_conf(2),
														typename insert_element<pga3::vector_conf(3), cl_null>::clist>::clist>::clist>::clist clist;

		typedef typename A::metric metric;
		typedef typename A::element_t element_t;

		explicit normalize_mv(const A &a_)
						: a(a_)
		{
		}

		auto operator()()
		{
			return ELEMENT(1)/sqrt((a & a).template element<0>()) * a;
		}
	private:
		A a;
	};

	// Even grade multivector norm
	template <class A>
	struct normalize_mv<A, 2> : public expression<normalize_mv<A>>
	{
		// The output configuration list is composed of even grades only
		typedef even_grade_clist clist;

		typedef typename A::metric metric;
		typedef typename A::element_t element_t;

		explicit normalize_mv(const A &a_)
						: a(a_)
		{
		}

		auto operator()()
		{
			auto& layout = pga3d<ELEMENT>::get();
			auto result = layout.even_grade_normalize(a);
			return result;
		}
	private:
		A a;
	};

	// Normalize other multivectors
	template <class A>
	struct normalize_mv<A, 0> : public expression<normalize_mv<A>>
	{
		typedef typename A::clist clist;
		typedef typename A::metric metric;
		typedef typename A::element_t element_t;

		explicit normalize_mv(const A &a_)
						: a(a_)
		{
		}

		auto operator()()
		{
			auto& layout = pga3d<ELEMENT>::get();
			return ELEMENT(1)/layout.norm(a) * a;
		}
	private:
		A a;
	};

	/// @brief The normalize of a PGA3 multivector
	/// @tparam E 	The expression class of the value
	/// @param x 		The expression
	/// @return 	The logarithm
	template <class E>
	auto normalize(const gaalet::expression<E>& x) const
	{
		return normalize_mv<E>(x)();
	}

	template <class E>
	auto norm(const gaalet::expression<E>& x) const
	{
		auto norm2 = (x*(~x)).template element<0>();
		if (::isclose(norm2, ELEMENT(0)))
		{
			// The element is ideal, so find the magnitude when joined with the origin, E0
			// (which is a representative Euclidean point) to compute the ideal norm
			// Note that e0 is the dual of E0
			auto x_joined_with_origin = dual((dual(x) ^ e0));
			norm2 = (x_joined_with_origin * (~x_joined_with_origin)).template element<0>();
		}
		return sqrt(norm2);
	}

	template <class E>
	auto even_grade_normalize(const gaalet::expression<E>& x) const
	{

		auto scalar_part = ::grade<0>(x).template element<0>();
		auto bivector_part = ::grade<2>(x);
		auto i_value = bivector_part.template element<pga3::I_CONF>();
		// Need to negate the j value since the value is stored for even permutation of the basis vectors,
		// but it is for an odd permutation.  This is an implementation detail.
		auto j_value = ELEMENT(-1) * bivector_part.template element<pga3::J_CONF>();
		auto k_value = bivector_part.template element<pga3::K_CONF>();
		auto dual_i_value = bivector_part.template element<pga3::DUAL_I_CONF>();
		auto dual_j_value = bivector_part.template element<pga3::DUAL_J_CONF>();
		auto dual_k_value = bivector_part.template element<pga3::DUAL_K_CONF>();
		auto pseudoscalar_part = ::grade<4>(x).template element<pseudoscalar_conf>();

		auto rotor_norm2 = scalar_part * scalar_part + i_value * i_value + j_value * j_value + k_value * k_value;
		if (::isclose(rotor_norm2, ELEMENT(0)))
		{
			// Ideal line
			auto norm2_ideal_line = dual_k_value * dual_k_value + dual_j_value * dual_j_value + dual_i_value * dual_i_value;
			auto inv_norm = ELEMENT(1) / sqrt(norm2_ideal_line);
			return make_motor(0, dual_k_value * inv_norm, dual_j_value * inv_norm, dual_i_value * inv_norm, 0, 0, 0, 0);
		}
		auto a = 1 / sqrt(rotor_norm2);
		auto b = (pseudoscalar_part * scalar_part -
		          (dual_k_value * k_value + dual_j_value * j_value + dual_i_value * i_value)) * a * a * a;
		return make_motor(a * scalar_part,
		             a * dual_k_value + b * k_value,
		             a * dual_j_value + b * j_value,
		             a * dual_i_value + b * i_value,
		             a * i_value,
		             a * j_value,
		             a * k_value,
		             a * pseudoscalar_part - b * scalar_part);

	}

	/// Square root of grade 0 and 2 multivectors only
	/// Objects to select the correct square root function at compile time using the grade
	// check for multivectors of that can be normalized
	// value=2 : even grade only
	// value=0 : scalar
	// value=-1 : not supported
	template<class A>
	struct sqrt_evaluation_type {
		static const int value =
						(check_even_grade<typename A::clist>::value) ? 2 :
						(check_scalar<typename A::clist>::value) ? 0 :
						-1;
	};


	/// Objects to select the correct logarithm function at compile time using the grade
	template<class A, int ET = sqrt_evaluation_type<A>::value>
	struct square_root : public expression<square_root<A>>
	{
		static_assert(ET!=-1, "no method for evaluating this type of multivector implemented");
	};

	/// PGA3 square root of even grade multivectors
	template <class A>
	struct square_root<A, 2> : public expression<square_root<A>>
	{
		typedef typename A::clist clist;
		typedef typename A::metric metric;
		typedef typename A::element_t element_t;

		explicit square_root(const A &a_)
						: a(a_)
		{
		}

		auto operator()()
		{
			auto& layout = pga3d<ELEMENT>::get();
			return layout.rotor_sqrt(a);
		}
	private:
		A a;
	};

	/// Scalar square root
	template <class A>
	struct square_root<A, 0> : public expression<square_root<A>>
	{
		// The output configuration list is composed of the scalar grade only
		typedef typename A::clist clist;
		typedef typename A::metric metric;
		typedef typename A::element_t element_t;

		explicit square_root(const A &a_)
						: a(a_)
		{
		}

		auto operator()()
		{
			return sqrt(a);
		}
	private:
		A a;
	};

	template <class E>
	auto rotor_sqrt(const gaalet::expression<E>& x) const
	{
		auto& layout = pga3d<ELEMENT>::get();
		return layout.normalize(one + x);
	}

	/// @brief The square root of a PGA3 multivector (only even grade multivectors and scalars supported)
	/// @tparam E 	The expression class of the value
	/// @param x 		The expression
	/// @return 	The square root
	template <class E>
	auto mv_sqrt(const gaalet::expression<E>& x) const
	{
		return square_root<E>(x)();
	}

	/// Inverse of a non-null blade
	template <class E>
	auto nonnull_blade_inverse(const E& v) const
	{
		auto reverse_v = ~v;
		assert((::grade<0>(reverse_v * v).template element<0>() != 0.));
		return 1./::grade<0>(reverse_v * v).template element<0>() * reverse_v;
	}

	///@brief Extract the coefficients of the point to a Euclidean point
	template <class E>
	std::tuple<ELEMENT> euc_point(const E& v) const
	{
		// Need to negate the x and z values since the values are stored for even permutation of the basis vectors,
		// but these are an odd permutation.  This is an implementation detail.
		return std::make_tuple(-v.template element<blade_conf(0, 3, 2)>(),
		                       v.template element<blade_conf(0, 1, 3)>(),
		                       -v.template element<blade_conf(0, 2, 1)>());
	}

	///@brief Extract the coefficients of the point to a Euclidean point
	template <class E>
	std::tuple<ELEMENT> down(const E& v) const
	{
		// Need to negate the x and z values since the values are stored for even permutation of the basis vectors,
		// but these are an odd permutation.  This is an implementation detail.
		return std::make_tuple(-v.template element<blade_conf(0, 3, 2)>(),
		                       v.template element<blade_conf(0, 1, 3)>(),
		                       -v.template element<blade_conf(0, 2, 1)>());
	}

	// Sandwich product
	template<class X, class V>
	auto sandwich(const gaalet::expression<X>& x, const gaalet::expression<V>& v) const
	{
		return v * x * ~v;
	}

	// Regressive product (implements a join in the common subspace)
	template<typename L, typename R>
	auto vee(const gaalet::expression<L> &l, const gaalet::expression<R> &r) const
	{
		return dual(dual(l) ^ dual(r));
	}

	template<typename L, typename R>
	inline
	auto line_from_points(const gaalet::expression<L> &start, const gaalet::expression<R> &end) const
	{
		return vee(start, end);
	}

	template<class T1, class T2, class T3>
	inline
	auto plane_from_points(const T1 &P1, const T2 &P2, const T3 &P3) const
	{
		return normalize(vee(vee(P1, P2), P3));
	}

	// Points are in a CCW direction
	template<class T1, class T2, class T3>
	inline
	auto normal_to_plane(const T1 &P1, const T2 &P2, const T3 &P3) const
	{
		return normalize(P1 & plane_from_points(P1, P2, P3));
	}

	// Commutator product
	template<typename L, typename R>
	auto commutator(const gaalet::expression<L> &l, const gaalet::expression<R> &r) const
	{
		return ELEMENT(0.5) * (l * r - r * l);
	}

	// Anti-commutator product
	template<typename L, typename R>
	auto anti_commutator(const gaalet::expression<L> &l, const gaalet::expression<R> &r) const
	{
		return ELEMENT(0.5) * (l * r + r * l);
	}

	template<class L, class R>
	auto meet(const gaalet::expression<L>& l, const gaalet::expression<R>& r) const
	{
		return l ^ r;
	}

	/// @brief The projection of the left object onto the right object's subspace
	/// @tparam L 	Class of the left object
	/// @tparam R 	Class of the right object
	/// @param l 		The left object
	/// @param r 		The right object
	/// @return  An object with the same grade as the left object
	template<class L, class R>
	auto project(const gaalet::expression<L>& l, const gaalet::expression<R>& r) const
	{
		return (l & r) * r;
	}

	/// @brief The rejection of the left object from the right object's subspace
	/// @tparam L 	Class of the left object
	/// @tparam R 	Class of the right object
	/// @param l 		The left object
	/// @param r 		The right object
	/// @return  An object with the same or a lesser grade than the right object
	template<class L, class R>
	auto reject(const gaalet::expression<L>& l, const gaalet::expression<R>& r) const
	{
		return (l & r);
	}

	/// @brief Reflect the left object in the given plane
	/// @tparam X 	Class of the left object
	/// @tparam P 	Class of the plane object
	/// @param x 		The object to reflect
	/// @param r 		The plane subspace of the object
	/// @return  An object with the same or a lesser grade than the right object
	template<class X, class P>
	auto reflect_in_plane(const gaalet::expression<X>& x, const gaalet::expression<P>& r) const
	{
		auto plane = ::part<PLANE_CONF>(r);
		return ELEMENT(-1) * plane * x * ~plane;
	}


	Rotor rotor_from_rotation_vector(const std::vector<ELEMENT> &rvec) const
	{
		ELEMENT r1 = rvec[0];
		ELEMENT r2 = rvec[1];
		ELEMENT r3 = rvec[2];

		ELEMENT angle_squared = r1 * r1 + r2 * r2 + r3 * r3;
		ELEMENT theta = sqrt(angle_squared);

		// Arbitrarily initialize the axis to be the z-axis in case the angle is 0
		Rotor axis = e12;
		if (theta > ELEMENT(0))
		{
			axis = r3 * e12 + r2 * e31 + r1 * e23;
			axis = ELEMENT(1)/ theta * axis;
		}
		else
		{
			theta = ELEMENT(0);   // Ensure that the derivatives are zero
		}

		typename pga3d<ELEMENT>::Scalar c = cos(ELEMENT(0.5) * theta);
		typename pga3d<ELEMENT>::Rotor s = sin(ELEMENT(0.5) * theta) * axis;
		typename pga3d<ELEMENT>::Rotor rotor = c - s;

		return rotor;
	}


	template <class E>
	std::vector<ELEMENT> rotation_vector_from_rotor(const gaalet::expression<E>& r) const
	{
		std::vector<ELEMENT> rvec = {ELEMENT(0), ELEMENT(0), ELEMENT(0)};
		if (::magnitude2(::grade<2>(r)).template element<0>() > 0 || ::grade<0>(r).template element<0>() != 1)
		{
			ELEMENT half_angle = acos(::grade<0>(r).template element<0>());
			typename pga3d<ELEMENT>::Rotor axis = eval(-1 * e12);    // Set an arbitrary axis if the rotation angle is 0
			if (half_angle != ELEMENT(0))
			{
				axis = eval(1. / sin(half_angle) * ::grade<2>(r));
			}
			Plane plane_of_rotation = ::grade<1>(axis * I3);
			rvec[0] = ELEMENT(2) * half_angle * plane_of_rotation.template element<pga3::vector_conf(1)>();
			rvec[1] = ELEMENT(2) * half_angle * plane_of_rotation.template element<pga3::vector_conf(2)>();
			rvec[2] = ELEMENT(2) * half_angle * plane_of_rotation.template element<pga3::vector_conf(3)>();
		}
		return rvec;
	}

/// @brief Compute the rotation rotor given an angle and an axis
/// @tparam E 	The type of the axis
/// @tparam ELEMENT 	The type of the element
/// @param axis 		The line around which the rotation is done
/// @param angle 		The angle of rotation
/// @return 	The rotation rotor
	template <class E>
	auto generate_rotation_rotor(const gaalet::expression<E>& axis, const ELEMENT angle) const
	{
		return one * cos(-angle * 0.5) + sin(-angle * 0.5) * normalize(axis);
	}

	template <class E>
	auto translator_from_direction_components(const std::vector<E>& direction) const
	{
		auto translator = one - ELEMENT(0.5) * (direction[0] * e01 + direction[1] * e02 + direction[2] * e03);
		return translator;
	}

	auto translator_from_direction(const Direction & direction) const
	{
		auto translator = one - ELEMENT(0.5) * (direction * I3);
		return translator;
	}

	/// @brief Compute the translation rotor given a Euclidean vector
/// @tparam ELEMENT 	The type of the element
/// @param t 		Direction (aka ideal point)
/// @return 	The translation rotor
	template <class E>
	Translator generate_translation_rotor(const gaalet::expression<E>& d) const
	{
		auto direction = ::part<DIRECTION_CONF>(d);
		auto translator = one + ELEMENT(-0.5) * direction * I3;
		return translator;
	}


	template<class L, class T>
	inline
	auto generate_translation_rotor(const L &line, const T &distance) const
	{
		return one + ELEMENT(0.5) * distance * normalize(line) * I;
	}

	Motor motor_from_rt_parameters(std::vector<ELEMENT> &extrinsic_parameters) const
	{
		std::vector<ELEMENT> direction_components({extrinsic_parameters[3], extrinsic_parameters[4], extrinsic_parameters[5]});
		typename pga3d<ELEMENT>::Translator translator = translator_from_direction_components<ELEMENT>(direction_components);
		typename pga3d<ELEMENT>::Rotor rotor = rotor_from_rotation_vector(extrinsic_parameters);

		return translator * rotor;
	}

	std::vector<ELEMENT> rt_parameters_from_motor(const Motor& motor) const
	{
		std::vector<ELEMENT> rt_parameters {ELEMENT(0), ELEMENT(0), ELEMENT(0), ELEMENT(0), ELEMENT(0), ELEMENT(0)};

		auto twist = log(motor);
		std::pair<ELEMENT, ELEMENT> study_number;
		typename pga3d<ELEMENT>::Bivector normalised_bivector;
		std::tie(study_number, normalised_bivector) = bivector_axis(twist);
		auto half_angle = study_number.first;

		// Separate the rotor and translator; the rotation part stays in log form since that is the Rodrigues vector too
		auto half_angle_axis = ::part<pga3::I_CONF, pga3::J_CONF, pga3::K_CONF>(half_angle * normalised_bivector);
		auto translator = ::part<0, pga3::DUAL_K_CONF, pga3::DUAL_J_CONF, pga3::DUAL_I_CONF>(motor * exp(ELEMENT(-1)* half_angle_axis));

		rt_parameters[0] = ELEMENT(-2) * half_angle_axis.template element<pga3::K_CONF>();  // x-axis
		rt_parameters[1] = ELEMENT(2) * half_angle_axis.template element<pga3::J_CONF>();   // y-axis; need to negate because value is stored negated
		rt_parameters[2] = ELEMENT(-2) * half_angle_axis.template element<pga3::I_CONF>();  // z-axis
		rt_parameters[3] = ELEMENT(-2) * translator.template element<pga3::DUAL_K_CONF>();
		rt_parameters[4] = ELEMENT(-2) * translator.template element<pga3::DUAL_J_CONF>();
		rt_parameters[5] = ELEMENT(-2) * translator.template element<pga3::DUAL_I_CONF>();

		return rt_parameters;
	}

	template <class E>
	std::vector<ELEMENT> rt_parameters_from_motor(const gaalet::expression<E>&& motor) const
	{
		return rt_parameters_from_motor(::part<MOTOR_CONF>(eval(motor)));
	}

	// Get a singleton object for the PGA3 space
	static const pga3d<ELEMENT>& get()
	{
		static const pga3d<ELEMENT>* pga = new pga3d<ELEMENT>();
		assert(pga != nullptr);
		return *pga;
	}


protected:
	/// Destructor should not be used since this is a singleton class
	~pga3d()
	{
		assert(false);
	}

}; // end class pga3d


#endif //PGA3D_H
