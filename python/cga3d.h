//
// Class to encapsulate operators and objects of 3d Conformal Geometric Algebra
//

#ifndef CGA3D_H
#define CGA3D_H

#pragma once

#include <array>
#include <cassert>
#include <cstdarg>
#include "gaalet.h"

// Convenience function and macros to generate the configuration list values
// Note that the order of the vectors in the blade_conf is not retained (because the result is just a bit string),
// but for consistency they will be shown in the order of the basis vectors in the blades
inline constexpr static gaalet::conf_t vector_conf(const int v = 0) { return v < 1 ? 0 : 1 << (v - 1); }

#define blade_conf(v1, v2, ...) (vector_conf(v1) | vector_conf(v2) | vector_conf( __VA_ARGS__ ))
#define blade4_conf(v1, v2, v3, v4) (vector_conf(v1) | vector_conf(v2) | vector_conf( v3 ) | vector_conf( v4 ))


template <typename ELEMENT>
class cga3 {

public:

	/// Define some configuration lists for common multivector types
#define EUC_VECTOR_CONF \
	vector_conf(1), vector_conf(2), vector_conf(3)
#define VECTOR_CONF \
	EUC_VECTOR_CONF, vector_conf(4), vector_conf(5)
#define ROTOR_CONF \
	0, blade_conf(1, 2), blade_conf(1, 3), blade_conf(2, 3)
#define TRANSLATOR_CONF \
	0, blade_conf(1, 4),  blade_conf(1, 5), blade_conf(2, 4), blade_conf(2, 5), blade_conf(3, 4), blade_conf(3, 5)
#define MOTOR_CONF \
	ROTOR_CONF, blade_conf(1, 4),  blade_conf(1, 5), blade_conf(2, 4), blade_conf(2, 5), \
	blade_conf(3, 4), blade_conf(3, 5), blade4_conf(1, 2, 3, 4), blade4_conf(1, 2, 3, 5)
#define LINE_CONF \
	blade_conf(1, 2, 4),  blade_conf(1, 2, 5), blade_conf(2, 3, 4), blade_conf(2, 3, 5), \
	blade_conf(1, 3, 4), blade_conf(1, 3, 5), blade_conf(1, 4, 5), blade_conf(2, 4, 5), blade_conf(3, 4, 5)
#define DUAL_LINE_CONF \
	blade_conf(1, 4), blade_conf(1, 5), blade_conf(2, 4), blade_conf(2, 5), blade_conf(3, 4), blade_conf(3, 5), \
	blade_conf(2, 3), blade_conf(1, 3), blade_conf(1, 2)
#define POINT_PAIR_CONF \
	DUAL_LINE_CONF, blade_conf(4, 5)
#define CIRCLE_CONF \
  LINE_CONF, blade_conf(1, 2, 3)
#define SPHERE_CONF \
  blade4_conf(1, 2, 3, 4), blade4_conf(1, 2, 3, 5), blade4_conf(2, 3, 4, 5), blade4_conf(2, 3, 4, 5), \
  blade4_conf(1, 3, 4, 5), blade4_conf(1, 3, 4, 5), blade4_conf(1, 2, 4, 5), blade4_conf(1, 2, 4, 5)
#define HOMOGENEOUS_POINT_CONF \
	blade_conf(1, 2, 4),  blade_conf(1, 2, 5), blade_conf(2, 3, 4), blade_conf(2, 3, 5), \
	blade_conf(1, 3, 4), blade_conf(1, 3, 5), blade_conf(1, 2, 3)

	///@brief Definition of the 5 dimensional conformal geometric algebra space
	typedef gaalet::algebra<gaalet::signature<4, 1>, ELEMENT> algebra;

	// The multivector type for the algebra
	template<gaalet::conf_t... elements>
	struct mv {
		typedef typename algebra::template mv<elements...>::type type;
	};

	/// Frequently used multivector types
	typedef typename mv<0x0>::type Scalar;
	typedef typename mv<EUC_VECTOR_CONF>::type EuclideanVector; // Used for Euclidean points
	typedef typename mv<VECTOR_CONF>::type Vector;  // Used for conformal points and dual spheres, and dual planes
	typedef typename mv<POINT_PAIR_CONF>::type PointPair;
	typedef typename mv<ROTOR_CONF>::type Rotor;
	typedef typename mv<TRANSLATOR_CONF>::type Translator;
	typedef typename mv<MOTOR_CONF>::type Motor;
	typedef typename mv<LINE_CONF>::type Line;
	typedef typename mv<DUAL_LINE_CONF>::type DualLine; // This corresponds to a PGA3 line
	typedef typename mv<CIRCLE_CONF>::type Circle;
	typedef typename mv<SPHERE_CONF>::type Sphere;
	typedef Sphere Plane;
	typedef Vector DualSphere;
	typedef Vector DualPlane; // This corresponds to a PGA3 plane
	typedef typename mv<HOMOGENEOUS_POINT_CONF>::type DualPoint;  // This corresponds to a PGA3 point

	/// The basis elements of the algebra

	/// Scalar
	Scalar one;

	/// Vectors
	typename mv<vector_conf(1)>::type e1;
	typename mv<vector_conf(2)>::type e2;
	typename mv<vector_conf(3)>::type e3;
	typename mv<vector_conf(4)>::type ep;
	typename mv<vector_conf(5)>::type em;

	/// Null vectors
	typename mv<vector_conf(4), vector_conf(5)>::type eo;
	typename mv<vector_conf(4), vector_conf(5)>::type einf;

	/// Bivectors used when representing motors
	typename mv<blade_conf(1, 2)>::type e12;    // rotation around z-axis
	typename mv<blade_conf(1, 3)>::type e31;    // rotation around y-axis
	typename mv<blade_conf(2, 3)>::type e23;    // rotation around x-axis
	typename mv<blade_conf(1, 4), blade_conf(1, 5)>::type e1inf; // translation in x-axis direction
	typename mv<blade_conf(2, 4), blade_conf(2, 5)>::type e2inf; // translation in y-axis direction
	typename mv<blade_conf(3, 4), blade_conf(3, 5)>::type e3inf; // translation in z-axis direction

	/// Other bivectors
	typename mv<blade_conf(1, 4), blade_conf(1, 5)>::type e1o;
	typename mv<blade_conf(2, 4), blade_conf(2, 5)>::type e2o;
	typename mv<blade_conf(3, 4), blade_conf(3, 5)>::type e3o;
	typename mv<blade_conf(4, 5)>::type E0;

	/// Pseudoscalar for Euclidean space
	typename mv<blade_conf(1, 2, 3)>::type I3;

	/// Other trivectors
	/// Lines are sums of: e23inf, e31inf, e12inf, e1oinf, e2oinf, e3oinf
	typename mv<blade_conf(2, 3, 4), blade_conf(2, 3, 5)>::type e23inf;
	typename mv<blade_conf(1, 3, 4), blade_conf(1, 3, 5)>::type e31inf;
	typename mv<blade_conf(1, 2, 4), blade_conf(1, 2, 5)>::type e12inf;
	typename mv<blade_conf(2, 3, 4), blade_conf(2, 3, 5)>::type e23o;
	typename mv<blade_conf(1, 3, 4), blade_conf(1, 3, 5)>::type e31o;
	typename mv<blade_conf(1, 2, 4), blade_conf(1, 2, 5)>::type e12o;
	typename mv<blade_conf(1, 4, 5)>::type e1oinf;
	typename mv<blade_conf(2, 4, 5)>::type e2oinf;
	typename mv<blade_conf(3, 4, 5)>::type e3oinf;

	/// Quadvectors
	typename mv<blade4_conf(1, 2, 3, 4), blade4_conf(1, 2, 3, 5)>::type e123inf;
	typename mv<blade4_conf(1, 2, 3, 4), blade4_conf(1, 2, 3, 5)>::type e123o;
	typename mv<blade4_conf(2, 3, 4, 5)>::type e23oinf;
	typename mv<blade4_conf(3, 1, 4, 5)>::type e31oinf;
	typename mv<blade4_conf(1, 2, 4, 5)>::type e12oinf;

	/// The pseudoscalar (representing all space)
	static constexpr gaalet::conf_t pseudoscalar_conf = {(1 << algebra::metric::dimension) - 1};
	typename mv<pseudoscalar_conf>::type I5;

	/// Constructor - initialises the basis elements
	cga3()
	{
		one = {ELEMENT(1.0)};
		e1 = {ELEMENT(1.0)};
		e2 = {ELEMENT(1.0)};
		e3 = {ELEMENT(1.0)};
		em = {ELEMENT(1.0)};
		ep = {ELEMENT(1.0)};
		eo = {eval(0.5 * (em - ep))};
		einf = {eval(em + ep)};
		e12 = {eval(cga3<ELEMENT>::e1 ^ e2)};
		e31 = {eval(cga3<ELEMENT>::e3 ^ e1)};
		e23 = {eval(cga3<ELEMENT>::e2 ^ e3)};
		e1inf = {eval(e1 ^ einf)};
		e2inf = {eval(e2 ^ einf)};
		e3inf = {eval(e3 ^ einf)};
		e1o = {eval(e1 ^ eo)};
		e2o = {eval(e2 ^ eo)};
		e3o = {eval(e3 ^ eo)};
		e23inf = {eval(e23 ^ einf)};
		e31inf = {eval(e31 ^ einf)};
		e12inf = {eval(e12 ^ einf)};
		e23o = {eval(e23 ^ eo)};
		e31o = {eval(e31 ^ eo)};
		e12o = {eval(e12 ^ eo)};
		e1oinf = {eval(e1 ^ eo ^ einf)};
		e2oinf = {eval(e2 ^ eo ^ einf)};
		e3oinf = {eval(e3 ^ eo ^ einf)};
		e123inf = {eval(e1 ^ e2 ^ e3 ^ einf)};
		e123o = {eval(e1 ^ e2 ^ e3 ^ eo)};
		e23oinf = {eval(e23 ^ eo ^ einf)};
		e31oinf = {eval(e31 ^ eo ^ einf)};
		e12oinf = {eval(e12 ^ eo ^ einf)};

		E0 = {eval(eo ^ einf)};
		I3 = {eval(e1 ^ e2 ^ e3)};
		I5 = {eval(e1 ^ e2 ^ e3 ^ ep ^ em)};

	}

	EuclideanVector coordinates_to_euc_point(const std::array<double, 3> & point) const
	{
		return point[0] * e1 + point[1] * e2 + point[2] * e3;
	}

	EuclideanVector coordinates_to_euc_point(const std::vector<double> & point) const
	{
		return point[0] * e1 + point[1] * e2 + point[2] * e3;
	}

	///@brief Raise a Euclidean point to a conformal point
	Vector up(const EuclideanVector & euc_vec) const
	{
		return eo + euc_vec + 0.5 * euc_vec * euc_vec * einf;
	}

	template <class E>
	static auto normalise(const E& v)
	{
		auto squared_mag = ::grade<0>((~v) * v).template element<0>();
		auto magnitude = sqrt(fabs(squared_mag));  // Need the fabs(), since the squared magnitude can be negative
		return 1./magnitude * v;
	}

	/// Inverse of a non-null blade
	template <class E>
	auto nonnull_blade_inverse(const E& v) const
	{
		auto reverse_v = ~v;
		assert((::grade<0>(reverse_v * v).template element<0>() != 0.));
		return 1./::grade<0>(reverse_v * v).template element<0>() * reverse_v;
	}

	/// Inverse homogenization - lowers a conformal point to a stereographic embedding
	Vector inv_homogenization(const Vector& vec) const
	{
		return vec * nonnull_blade_inverse((-1. * vec) & einf);
	}

	///@brief Lower a conformal point to a Euclidean point
	EuclideanVector down(const Vector& conformal_vec) const
	{
		// Rejection of the inverse homogenized vector by E0 (note that 1/E0 == E0, so * is used instead of /)
		return (inv_homogenization(conformal_vec) ^ E0) * E0;
	}

	template<class A>
	auto dual(const gaalet::expression<A>& a) const
	{ // gaalet::dual has a bug for CGA3 (it gets the sign wrong), so implement it as a regular function
		return a & (-1. * I5);
	}

	template<class L, class R>
	auto meet(const gaalet::expression<L>& l, const gaalet::expression<R>& r) const
	{
		return I5 * ((I5 * l) ^ (I5 * r));
	}

	/// @brief Project a multivector l on the sub-space defined by the multivector r
	/// @param	l	The multivector to project
	/// @param	r		The multivector defining the sub-space
	/// @return The scalar value
	template<class L, class R>
	auto project(const gaalet::expression<L>& l, const gaalet::expression<R>& r)
	{
		return (l & r) * inverse(r);
	}

	/// Normalises a conformal point so that it has an inner product of -1 with einf
	template<class E>
	auto normalise_n_minus_1(const gaalet::expression<E>& mv) const
	{
		auto scale = ::grade<0>(mv & einf).template element<0>();
		if (scale != 0)
			return (-1 / scale) * mv;

		assert(false && "Couldn't normalise the conformal point");
		return 1 * mv;    // Multiply by 1 to satisfy the compiler's auto deduction mechanism
	}

	/// Extracts the end points of a point pair bivector
	/// translated from clifford.tools.g3c.point_pair_to_end_points()
	template<class E>
	auto point_pair_to_end_points(const gaalet::expression<E>& pp) const
	{
		auto beta = sqrt(fabs((pp * pp).template element<0>()));
		auto F = ELEMENT(1.)/ beta * pp;
		auto P = ELEMENT(0.5) * F + ELEMENT(0.5) * one;
		auto P_twiddle = ELEMENT(-0.5) * F + ELEMENT(0.5) * one;
		auto A = normalise_n_minus_1(ELEMENT(-1.) * (P_twiddle * (pp & einf)));
		auto B = normalise_n_minus_1((P * (pp & einf)));

		return std::make_pair(A, B);
	}


	// Get a singleton object for the CGA3 space
	static const cga3<ELEMENT>& get()
	{
		static cga3<ELEMENT>* cga = nullptr;

		if (cga == nullptr)
		{
			// Allocate a singleton object for the CGA3 space
			cga = new cga3<ELEMENT>();
		}
		return *cga;
	}


protected:
	/// Destructor should not be used since this is a singleton class
	~cga3()
	{
		assert(false);
	}

}; // end class cga3


#endif //CGA3D_H
