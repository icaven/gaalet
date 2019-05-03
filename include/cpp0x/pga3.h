#include "gaalet.h"
#include <cmath>
#include <iostream>

#pragma once

using namespace gaalet;

// Projective Geometric Algebra 3d
// For a description of these operations, see:
// Gunn, Charles, “Projective geometric algebra: A new framework for doing euclidean geometry”,
// https://export.arxiv.org/abs/1901.05873
//
// Also, for a javascript implementation and description, see:
// https://github.com/enkimute/ganja.js

namespace pga3
{

struct space {
    typedef gaalet::algebra<signature<3, 0, 1>, double> algebra;

    template <gaalet::conf_t... elements> struct mv {
        typedef typename algebra::mv<elements...>::type type;
    };
};

// Convenience function and macro to generate the configuration list values
// Note that the null basis e0 has a bit representation that is expected by the gaalet implementation to square to 0
inline constexpr conf_t vector_conf(const int v=-1) { return v < 0 ? 0 : 1 << ((v - 1) < 0 ? 3 : (v - 1)); };
#define blade_conf(v1, v2, ...) (pga3::vector_conf(v1) | pga3::vector_conf(v2) | pga3::vector_conf( __VA_ARGS__ ) )


// The basis elements of the dual algebra

// Scalar
const space::mv<0x0>::type one = { 1.0 };

// Vectors (representing planes)
const space::mv<vector_conf(0)>::type e0 = { 1.0 };
const space::mv<vector_conf(1)>::type e1 = { 1.0 };
const space::mv<vector_conf(2)>::type e2 = { 1.0 };
const space::mv<vector_conf(3)>::type e3 = { 1.0 };

// Bivectors (representing lines)
const space::mv<blade_conf(0, 1)>::type e01 = (e0 ^ e1);
const space::mv<blade_conf(0, 2)>::type e02 = (e0 ^ e2);
const space::mv<blade_conf(0, 3)>::type e03 = (e0 ^ e3);
const space::mv<blade_conf(1, 2)>::type e12 = (e1 ^ e2);
const space::mv<blade_conf(1, 3)>::type e13 = (e1 ^ e3);
const space::mv<blade_conf(2, 3)>::type e23 = (e2 ^ e3);

// Alternative names for biquaternion components
constexpr conf_t DUAL_I_CONF = blade_conf(0, 3);
constexpr conf_t DUAL_J_CONF = blade_conf(0, 2);
constexpr conf_t DUAL_K_CONF = blade_conf(0, 1);
constexpr conf_t I_CONF = blade_conf(1, 2);
constexpr conf_t J_CONF = blade_conf(1, 3);
constexpr conf_t K_CONF = blade_conf(2, 3);
const space::mv<I_CONF>::type i = e12;
const space::mv<J_CONF>::type j = -1 * e13;
const space::mv<K_CONF>::type k = e23;
const space::mv<DUAL_K_CONF>::type dual_k = e01;
const space::mv<DUAL_J_CONF>::type dual_j = e02;
const space::mv<DUAL_I_CONF>::type dual_i = e03;


// The tri-vectors (representing points)
constexpr conf_t E0_CONF = blade_conf(1, 2, 3);
constexpr conf_t EX_CONF = blade_conf(0, 1, 2);
constexpr conf_t EY_CONF = blade_conf(0, 1, 3);
constexpr conf_t EZ_CONF = blade_conf(0, 2, 3);

const space::mv<blade_conf(1, 2, 3)>::type e123 = (e1 ^ e2 ^ e3);
const space::mv<blade_conf(0, 1, 2)>::type e012 = (e0 ^ e1 ^ e2);
const space::mv<blade_conf(0, 1, 3)>::type e013 = (e0 ^ e1 ^ e3);
const space::mv<blade_conf(0, 2, 3)>::type e023 = (e0 ^ e2 ^ e3);
// Alternative names
const space::mv<blade_conf(1, 2, 3)>::type E0 = e123; // This is also the pseudoscalar for the Euclidean subspace
const space::mv<blade_conf(0, 1, 2)>::type E3 = e012;
const space::mv<blade_conf(0, 1, 3)>::type E2 = e013;
const space::mv<blade_conf(0, 2, 3)>::type E1 = e023;
const space::mv<blade_conf(0, 1, 2)>::type EX = e012;
const space::mv<blade_conf(0, 1, 3)>::type EY = e013;
const space::mv<blade_conf(0, 2, 3)>::type EZ = e023;

// The pseudoscalar (representing all space), with the null vector at the end to match the order
// of the configuration bits
const gaalet::conf_t pseudoscalar_conf = (1 << space::algebra::metric::dimension)-1;
const space::mv<pseudoscalar_conf>::type I = (e1 ^ e2 ^ e3 ^ e0);

// Types for common geometric entities
typedef space::mv<E0_CONF, EX_CONF, EY_CONF, EZ_CONF>::type Point_t;
typedef space::mv<DUAL_K_CONF, DUAL_J_CONF, DUAL_I_CONF, I_CONF, J_CONF, K_CONF>::type Line_t;
typedef space::mv<vector_conf(1), vector_conf(2), vector_conf(3), vector_conf(0)>::type Plane_t;

//
// Specialized version of the dual function for the PGA algebra 
// TODO: modify the library dual function to support this algebra, since this version
// is only a slight modification of the library function
// Implements the Poincare duality
// See: section 2.3.1.3 of "Geometry, Kinematics, and Rigid Body Mechanics in Cayley-Klein Geometries"
// The sign change is needed for some blades since the basis vectors in dual configuration 
// are ordered to be in the canonical order instead of being permuted (as described by the referenced section)

namespace detail
{
template<conf_t I, typename list, typename colist = cl_null>
struct dual_list
{
   typedef typename dual_list<I, typename list::tail, typename insert_element< I ^ list::head, colist>::clist>::clist clist;
};
template<conf_t I, typename colist>
struct dual_list<I, cl_null, colist>
{
   typedef colist clist;
};

template<class A>
struct dual : public expression<dual<A>>
{
   static const conf_t I = Power<2, A::metric::dimension>::value-1;

   typedef typename dual_list<I, typename A::clist>::clist clist;

   typedef typename A::metric metric;
   
   typedef typename A::element_t element_t;

   dual(const A& a_)
      :  a(a_)
   { }

   template<conf_t conf>
   element_t element() const {
      return (search_element<conf, clist>::index>=clist::size) ? 0.0 : a.template element< I ^ conf >()
                                                                     * ((conf == blade_conf(0, 2) || conf == blade_conf(1, 3) || 
                                                                         conf == vector_conf(2)  || conf == blade_conf(0, 1, 3) ||
                                                                         conf == vector_conf(0) || conf == blade_conf(1, 2, 3)
                                                                         ? -1 : 1)
                                                                     );
   }

protected:
   const A a;
};

} // end namespace detail

/// Dual of a multivector.
/**  Specialized version for PGA3
  */
/// \ingroup ga_ops

#undef dual
template <class A> inline
detail::dual<A>
dual(const gaalet::expression<A>& a) {
   return detail::dual<A>(a);
}


//
// Functions to generate common entities
//

// Eulidean point as a homogeneous point
template <typename T> auto make_point(T x, T y, T z)
{
    return E0 + x * EX + y * EY + z * EZ;
}

inline space::algebra::element_t Point_x(const Point_t& p) {
    return p.template element<EX_CONF>();
};
inline space::algebra::element_t Point_y(const Point_t& p) {
    return p.template element<EY_CONF>();
};
inline space::algebra::element_t Point_z(const Point_t& p) {
    return p.template element<EZ_CONF>();
};

// Ideal point
template <typename T> auto make_ideal_point(T x, T y, T z)
{
    return space::algebra::element_t(0) * E0 + x * EX + y * EY + z * EZ;
}


// Lines can be defined by Plücker coordinates
template <typename TX, typename TY, typename TZ, typename DX, typename DY, typename DZ>
auto make_line(TX px, TY py, TZ pz, DX dx, DY dy, DZ dz)
{
    return normalize(px * dual_k + py * dual_j + pz * dual_i + dx * i + dy * j + dz * k);
}


// TODO: Uncomment and test this function 
//template <typename CL>
//auto make_line(const Point_t& origin, const space::mv<CL>& direction)
//{
//    auto extracted_direction = ::grade<2>(direction);
//    return Point_x(origin) * dual_k + Point_y(origin) * dual_j + Point_z(origin) * dual_i 
//            + extracted_direction.template element<I_CONF>() * i
//            + extracted_direction.template element<J_CONF>() * j 
//            + extracted_direction.template element<K_CONF>() * k;
//}


// Four values define a plane
template <typename T> auto make_plane(T a, T b, T c, T d)
{
    return normalize(d * e0 + a * e1 + b * e2 + c * e3);
}

//
// Geometric operations
//

template <typename X> inline auto polar(const gaalet::expression<X>& x)
{
//    return x * I;
    // Need to add the zero scalar so that the configuration list won't be empty after some geometric products
    return 0.0 * one + x * I;
    // The next line may be removed after the current strategy of adding the zero scalar is shown to work
//    return -(I * x);  // Wanted to do: (x * I) but ran into a template expansion ambiguity compiler error, so use anti-symmetric property
}

// Versor sandwich product
template<typename L, typename R> inline
auto sandwich(const gaalet::expression<L>& l, const gaalet::expression<R>& r)
{
   return (r*l*(~r));
}

// Joins
template<typename L, typename R> inline
auto vee(const gaalet::expression<L>& l, const gaalet::expression<R>& r)
{
    return pga3::dual(pga3::dual(l) ^ pga3::dual(r));
}

template<typename L, typename R> inline
auto line_from_points(const gaalet::expression<L>& start, const gaalet::expression<R>& end)
{
    return vee(start, end);
}

template <class T1, class T2, class T3> auto plane_from_points(const T1& P1, const T2& P2, const T3& P3)
{
    return normalize(vee(vee(P1, P2), P3));
}

template <class TP, class TL> auto plane_from_point_and_line(const TP& P, const TL& L)
{
    return normalize(vee(P, L));
}

// Meets
template <class T1, class T2> auto line_from_planes(const T1& p1, const T2& p2)
{
    return normalize(p1) ^ normalize(p2);
}
template <class T1, class T2, class T3> auto point_from_planes(const T1& p1, const T2& p2, const T3& p3)
{
    return normalize(p1 ^ p2 ^ p3);
}

template <class TP, class TL> auto point_from_line_and_plane(const TL& L, const TP& P)
{
    return normalize(L ^ P);
}

// Points are in a CCW direction
template <class T1, class T2, class T3> auto normal_to_plane(const T1& P1, const T2& P2, const T3& P3) {
    return -1. * normalize(P1 & plane_from_points(P1, P2, P3));
}

//
// Rigid motions
//

// Rotation around an axis
template <class A, class T> auto rotor(const A& axis, const T& angle)
{
    return one * cos(-angle * 0.5) + sin(-angle * 0.5) * normalize(axis);
};

// Translation in the direction of a line
template <class L, class T> auto translator(const L& line, const T& distance)
{
    return one + 0.5 * distance * normalize(line) * I;
}

// Motor operation implements a rotation followed by a translation, which is similar to
// a screw operation, which is a rotation around a line while travelling a distance along the line
template <class L, class T, typename P = space::algebra::element_t>
auto motor(const L& line, const T& distance, const P angle)
{
    return translator(line, distance) * rotor(line, angle);
}

//template <class L, class T, typename P = space::algebra::element_t>
//auto screw(const L& line, const T& distance, const P pitch)
//{
//    return exp(0.5 * distance * (normalize(line) + pitch * normalize(line) * I));
//}

// Alternative form of the screw operation, to be used for comparison
//template <class L, class T, typename P = space::algebra::element_t>
//auto screw2(const L& line, const T& distance, const P pitch)
//{
//    return exp(0.5 * distance * (one + pitch * I) * normalize(line));
//}

} // end of namespace pga3
