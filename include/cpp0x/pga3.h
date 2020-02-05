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
// Note that the order of the vectors in the blade_conf is not retained (because the result is just a bit string),
// but for consistency they will be shown in the order of the basis vectors in the blades
inline constexpr conf_t vector_conf(const int v=-1) { return v < 0 ? 0 : 1 << v; };
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
const space::mv<blade_conf(1, 3)>::type e31 = (e3 ^ e1);
const space::mv<blade_conf(2, 3)>::type e23 = (e2 ^ e3);

// Alternative names for biquaternion components
constexpr conf_t DUAL_I_CONF = blade_conf(0, 3);
constexpr conf_t DUAL_J_CONF = blade_conf(0, 2);
constexpr conf_t DUAL_K_CONF = blade_conf(0, 1);
constexpr conf_t I_CONF = blade_conf(1, 2);
constexpr conf_t J_CONF = blade_conf(3, 1);
constexpr conf_t K_CONF = blade_conf(2, 3);
const space::mv<I_CONF>::type i = e12;      // z-axis      - see Section 7.4.3 of [Gunn2011]
const space::mv<J_CONF>::type j = e31;      // y-axis
const space::mv<K_CONF>::type k = e23;      // x-axis
const space::mv<DUAL_K_CONF>::type dual_k = e01;
const space::mv<DUAL_J_CONF>::type dual_j = e02;
const space::mv<DUAL_I_CONF>::type dual_i = e03;


// The tri-vectors (representing points)
constexpr conf_t E0_CONF = blade_conf(1, 2, 3);
constexpr conf_t E1_CONF = blade_conf(0, 3, 2);
constexpr conf_t E2_CONF = blade_conf(0, 1, 3);
constexpr conf_t E3_CONF = blade_conf(0, 2, 1);


const space::mv<E0_CONF>::type e123 = (e1 ^ e2 ^ e3);
const space::mv<E1_CONF>::type e032 = (e0 ^ e3 ^ e2);
const space::mv<E2_CONF>::type e013 = (e0 ^ e1 ^ e3);
const space::mv<E3_CONF>::type e021 = (e0 ^ e2 ^ e1);

    // Alternative names
const space::mv<E0_CONF>::type E0 = e123; // This is also the pseudoscalar for the Euclidean subspace
const space::mv<E1_CONF>::type E1 = e032;
const space::mv<E2_CONF>::type E2 = e013;
const space::mv<E3_CONF>::type E3 = e021;

// The pseudoscalar (representing all space), with the null vector at the end to match the order
// of the configuration bits
const gaalet::conf_t pseudoscalar_conf = (1 << space::algebra::metric::dimension)-1;
const space::mv<pseudoscalar_conf>::type I = (e0 ^ e1 ^ e2 ^ e3);

// Types for common geometric entities
typedef space::mv<E0_CONF, E1_CONF, E2_CONF, E3_CONF>::type Point_t;
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
template<conf_t I_conf, typename list, typename colist = cl_null>
struct dual_list
{
   typedef typename dual_list<I_conf, typename list::tail, typename insert_element< I_conf ^ list::head, colist>::clist>::clist clist;
};
template<conf_t I_conf, typename colist>
struct dual_list<I_conf, cl_null, colist>
{
   typedef colist clist;
};

template<class A>
struct dual : public expression <detail::dual<A >>
{
   static const conf_t I_conf = Power<2, A::metric::dimension>::value-1;

   typedef typename dual_list<I_conf, typename A::clist>::clist clist;

   typedef typename A::metric metric;
   
   typedef typename A::element_t element_t;

   dual(const A& a_)
      :  a(a_)
   { }

   template<conf_t conf>
   element_t element() const {
       // The odd permutation blades are stored negated, so change the sign when computing the dual
      return (search_element<conf, clist>::index>=clist::size) ? 0.0 : a.template element< I_conf ^ conf >()
              //                                                                even                            odd
                                                                     * ((conf == blade_conf(0, 2) || conf == blade_conf(3, 1) ||
                                                                         conf == vector_conf(1)  || conf == blade_conf(0, 3, 2) ||
                                                                         conf == vector_conf(3) || conf == blade_conf(0, 2, 1)
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

// Euclidean point as a homogeneous point
template <typename X, typename Y, typename Z> inline
auto make_point(X x, Y y, Z z)
{
    return E0 + static_cast<pga3::space::algebra::element_t>(x) * E1 +
                static_cast<pga3::space::algebra::element_t>(y) * E2 +
                static_cast<pga3::space::algebra::element_t>(z) * E3;
}

// Ideal point
template <typename X, typename Y, typename Z> inline
auto make_ideal_point(X x, Y y, Z z)
{
    return static_cast<pga3::space::algebra::element_t>(x) * E1 +
           static_cast<pga3::space::algebra::element_t>(y) * E2 +
           static_cast<pga3::space::algebra::element_t>(z) * E3;
}

template<class A>
struct Point : public gaalet::expression<Point<A>>
{
    typedef typename A::clist clist;

    typedef typename A::metric metric;

    typedef typename A::element_t element_t;

    explicit Point(const A& a_)
            :  a(a_)
    { }

    template<conf_t conf>
    element_t element() const {
        return a.template element<conf>();
    }

    // Need to negate the x and z values since the value is stored for even permutation of the basis vectors,
    // but they are for an odd permutation.  This is an implementation detail.

    element_t origin() const {
        return a.template element<E0_CONF>();
    }

    // x direction
    element_t x() const {
        return -a.template element<E1_CONF>();
    }

    // y direction
    element_t y() const {
        return a.template element<E2_CONF>();
    }

    // z direction
    element_t z() const {
        return -a.template element<E3_CONF>();
    }

    Point_t normalized() const {
        return origin() * E0 + x() * E1 + y() * E2 + z() * E3;
    }


protected:
    const A a;
};

// Ensure that the point is normalized
template <class A> inline
Point<gaalet::unit<A>>
as_point(const gaalet::expression<A>& a) {
    gaalet::unit<A> na = ::normalize(a);
    return Point<decltype(na)>(na);
}

template <class A>
std::ostream& operator << (std::ostream& out, const Point<A>& p)
{
    return out << "(x: " << p.x() << ", y: " << p.y() << ", z: " << p.z() << ", O: " << p.origin() << ")";
}

// Need to negate the x and z values since the value is stored for even permutation of the basis vectors,
// but they are for an odd permutation.  This is an implementation detail.
inline space::algebra::element_t Point_x(const Point_t& p) {
    return -p.template element<E1_CONF>();
};
inline space::algebra::element_t Point_y(const Point_t& p) {
    return p.template element<E2_CONF>();
};
inline space::algebra::element_t Point_z(const Point_t& p) {
    return -p.template element<E3_CONF>();
};


// Lines can be defined by Plücker coordinates
template <typename TX, typename TY, typename TZ, typename DX, typename DY, typename DZ> inline
auto make_line(TX px, TY py, TZ pz, DX dx, DY dy, DZ dz)
{
    return normalize(dx * dual_k + dy * dual_j + dz * dual_i + px * k + py * j + pz * i);
}

template<class A>
struct Line : public gaalet::expression<Line<A>>
{
    typedef typename A::clist clist;

    typedef typename A::metric metric;

    typedef typename A::element_t element_t;

    explicit Line(const A& a_)
            :  a(a_)
    { }

    template<conf_t conf>
    element_t element() const {
        return a.template element<conf>();
    }

    // Need to negate the j value since the value is stored for even permutation of the basis vectors,
    // but it is for an odd permutation.  This is an implementation detail.

    element_t i() const {
        return a.template element<I_CONF>();
    }

    element_t j() const {
        return -a.template element<J_CONF>();
    }

    element_t k() const {
        return a.template element<K_CONF>();
    }

    element_t dual_i() const {
        return a.template element<DUAL_I_CONF>();
    }

    element_t dual_j() const {
        return a.template element<DUAL_J_CONF>();
    }

    element_t dual_k() const {
        return a.template element<DUAL_K_CONF>();
    }

    Line_t normalized() const {
        return make_line(k(), j(), i(), dual_k(), dual_j(), dual_i());
    }

protected:
    const A a;
};

template <class A>
std::ostream& operator << (std::ostream& out, const Line<A>& l)
{
    return out << "(i: " << l.i()
               << ", j: " << l.j()
               << ", k: " << l.k()
               << ", di: " << l.dual_i()
               << ", dj: " << l.dual_j()
               << ", dk: " << l.dual_k()
               << ")";
}

// Ensure that the line is normalized
template <class A> inline
Line<gaalet::unit<A>>
as_line(const gaalet::expression<A>& a) {
    gaalet::unit<A> na = ::normalize(a);
    return Line<decltype(na)>(na);
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

namespace pga3 {
// These functions may be used to print out PGA3 multivectors in a more descriptive way than just
// using the functions in streaming.h.
// @todo: Find a way to call these functions for gaalet::expressions that contain pga3 multivectors


// These are specific to PGA3, and are ordered so that the configuration bit pattern can index into it
static std::string basis_vector_names[16] = {"    ", "e0  ", "e1  ", "e01 ", "e2  ", "e02 ", "e12 ", "e021 ", "e3  ",
                                             "e03 ", "e31 ", "e013 ", "e23 ", "e032 ", "e123 ", "e0123"};

// PGA3 multivector streaming

    template <typename clist>
    struct UnpackElementsToStream {
        template<class E, class T>
        static void unpack(std::basic_ostream<E, T> &os, const gaalet::multivector<clist,
                pga3::space::algebra::metric, pga3::space::algebra::element_t>& e,
                bool previous_non_zero=false) {
            bool found_non_zero = previous_non_zero;
            if (e.template element<clist::head>() != 0) {
                found_non_zero = true;
                if (clist::head == 0) {
                    os << std::right << std::setw(8) << e.template element<0>();
                }
                else if (clist::head == 7 || clist::head == 10 || clist::head == 13) {
                    // The odd permutation blades are stored with negated values, so change the sign
                    os << std::right << std::setw(8) << -1 * e.template element<clist::head>() << basis_vector_names[clist::head];
                }
                else {
                    os << std::right << std::setw(8) << e.template element<clist::head>() << basis_vector_names[clist::head];
                }
            }
            UnpackElementsToStream<typename clist::tail>::unpack(os, e, found_non_zero);
        }

    };

    template<>
    struct UnpackElementsToStream<gaalet::cl_null> {
        template<class E, class T>
        static void unpack(std::basic_ostream<E, T> &os, const gaalet::multivector<gaalet::cl_null,
                pga3::space::algebra::metric, pga3::space::algebra::element_t>& , bool previous_non_zero=false)
        {
            if (!previous_non_zero) {
                os << std::right << std::setw(8) << "0";
            }
        }
    };

} // end of namespace pga3

namespace std {

    template<typename CL>
    std::ostream &operator<<(std::ostream &os, const gaalet::multivector<CL,
            pga3::space::algebra::metric, pga3::space::algebra::element_t> &e) {
        pga3::UnpackElementsToStream<CL>::unpack(os, e);
        return os;
    }

    template<class E, class T, class CL>
    std::basic_ostream<E, T> &operator<<(std::basic_ostream<E, T> &&os,
                                         const gaalet::expression<gaalet::multivector<CL, pga3::space::algebra::metric,
                                                 pga3::space::algebra::element_t>> &e) {
        pga3::UnpackElementsToStream<CL>::unpack(os, e);
        return os;
    }

} // end of namespace std
