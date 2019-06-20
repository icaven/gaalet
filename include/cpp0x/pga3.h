#pragma once


#include "gaalet.h"
#include <cmath>
#include <cfloat>
#include <iostream>


using namespace gaalet;

// Projective Geometric Algebra 3d
// For a description of these operations, see [Gunn2019]:
// Gunn, Charles, “Projective geometric algebra: A new framework for doing euclidean geometry”,
// https://export.arxiv.org/abs/1901.05873
//
// and: 
// Gunn, Charles, "Geometry, Kinematics, and Rigid Body Mechanics in Cayley-Klein Geometries" [Gunn2011]


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
const space::mv<I_CONF>::type i = e12;      // z-axis      - see Section 7.4.3 of [Gunn2011]
const space::mv<J_CONF>::type j = -1 * e13; // y-axis
const space::mv<K_CONF>::type k = e23;      // x-axis
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
typedef space::mv<0, DUAL_K_CONF, DUAL_J_CONF, DUAL_I_CONF, I_CONF, J_CONF, K_CONF, pseudoscalar_conf>::type Motor_t;

//
// Specialized version of the dual function for the PGA algebra 
// TODO: modify the library dual function to support this algebra, since this version
// is only a slight modification of the library function
// Implements the Poincare duality
// See: section 2.3.1.3 of Gunn, Charles, "Geometry, Kinematics, and Rigid Body Mechanics in Cayley-Klein Geometries" [Gunn2011]
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
struct dual : public expression <detail::dual<A >>
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
template <typename T> inline 
auto make_point(T x, T y, T z)
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
template <typename T> inline 
auto make_ideal_point(T x, T y, T z)
{
    return space::algebra::element_t(0) * E0 + x * EX + y * EY + z * EZ;
}

// Lines can be defined by Plücker coordinates
template <typename TX, typename TY, typename TZ, typename DX, typename DY, typename DZ> inline
auto make_line(TX px, TY py, TZ pz, DX dx, DY dy, DZ dz)
{
    return normalize(px * dual_k + py * dual_j + pz * dual_i + dx * i + dy * j + dz * k);
}

// Four values define a plane
template <typename T> inline 
auto make_plane(T a, T b, T c, T d)
{
    return normalize(d * e0 + a * e1 + b * e2 + c * e3);
}

//
// Geometric operations
//

template <typename X> inline 
auto polar(const gaalet::expression<X>& x)
{
//    return x * I;
    // Need to add the zero scalar so that the configuration list won't be empty after some geometric products
    return 0.0 * one + x * I;
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

template <class T1, class T2, class T3> inline
auto plane_from_points(const T1& P1, const T2& P2, const T3& P3)
{
    return normalize(vee(vee(P1, P2), P3));
}

template <class TP, class TL> inline
auto plane_from_point_and_line(const TP& P, const TL& L)
{
    return normalize(vee(P, L));
}

// Meets
template <class T1, class T2> inline
auto line_from_planes(const T1& p1, const T2& p2)
{
    return normalize(p1) ^ normalize(p2);
}
template <class T1, class T2, class T3> inline
auto point_from_planes(const T1& p1, const T2& p2, const T3& p3)
{
    return normalize(p1 ^ p2 ^ p3);
}

template <class TP, class TL> inline
auto point_from_line_and_plane(const TL& L, const TP& P)
{
    return normalize(L ^ P);
}

// Points are in a CCW direction
template <class T1, class T2, class T3> inline
auto normal_to_plane(const T1& P1, const T2& P2, const T3& P3) {
    return -1. * normalize(P1 & plane_from_points(P1, P2, P3));
}

//
// Rigid motions
//

// Rotation around an axis
template <class A, class T> inline
auto rotor(const A& axis, const T& angle)
{
    return one * cos(-angle * 0.5) + sin(-angle * 0.5) * normalize(axis);
};

// Translation in the direction of a line
template <class L, class T> inline
auto translator(const L& line, const T& distance)
{
    return one + 0.5 * distance * normalize(line) * I;
}

// Motor operation implements a rotation followed by a translation, which is similar to
// a screw operation, which is a rotation around a line while travelling a distance along the line
template <class L, class T, typename P = space::algebra::element_t> inline
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
template <class L, class T, typename P = space::algebra::element_t>
auto screw(const L& line, const T& distance, const P pitch)
{
    return 0.5 * distance * (one + pitch * I) * normalize(line);
}

inline bool isclose(double a, double b, double rtol=1e-05, double atol=DBL_EPSILON) {
    return fabs(a - b) <= (atol + rtol * fabs(b));
}

/// Return the Study number and the axis (a normalized bivector) associated with the given bivector
// See: Section 7.7 of Gunn, Charles, "Geometry, Kinematics, and Rigid Body Mechanics in Cayley-Klein Geometries" [Gunn2011]
// Normalize the bivector by multiplying it by the inverse of the square root of the associated Study number z, and the result
// is B, which is the axis of the given bivector
std::pair<std::pair<double, double>, Line_t> inline 
bivector_axis(pga3::Line_t bivector) {
    auto z = eval(bivector * (~bivector));
    double a = ::grade<0>(z).template element<0>();
    auto axis = eval(bivector);   // In the case that the inverse sqrt of z doesn't exist (when a == 0), the B is equal to the bivector
    double b = 0.0;
    std::pair<double, double> study_number = {0.0, 0.0};
    if (!isclose(a, 0.0)) {
        auto sqrt_a = sqrt(a);
        b = ::grade<4>(z).template element<pseudoscalar_conf>();
        auto inv_sqrt_z = sqrt_a/a * one - ((b)/(2 * a * sqrt_a))*I;
        axis = inv_sqrt_z * axis;
        study_number = std::make_pair(sqrt_a, b/(2.0 * sqrt_a));
        return std::make_pair(study_number, axis);
    }
    return std::make_pair(study_number, axis);
}

} // end of namespace pga3
