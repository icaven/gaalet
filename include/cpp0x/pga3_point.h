#ifndef __GAALET_PGA3_POINT_H
#define __GAALET_PGA3_POINT_H

#pragma once

#include "grade_check.h"
#include "pga3.h"
#include "pga3_normalize.h"

//
// PGA3 Point support
//

namespace pga3 {

    /// @brief Functional interface to making a Euclidean point as a normalized homogeneous point
    template<typename X, typename Y, typename Z>
    inline
    auto make_point(const X x, const Y y, const Z z) {
        return E0 + static_cast<space::algebra::element_t>(x) * E1 +
               static_cast<space::algebra::element_t>(y) * E2 +
               static_cast<space::algebra::element_t>(z) * E3;
    }

    /// @brief Functional interface to making an ideal point (non-Euclidean) - not normalized
    template<typename X, typename Y, typename Z>
    inline
    auto make_ideal_point(const X x, const Y y, const Z z) {
        return static_cast<space::algebra::element_t>(x) * E1 +
               static_cast<space::algebra::element_t>(y) * E2 +
               static_cast<space::algebra::element_t>(z) * E3;
    }

    /// @brief Accessors for the components of a point
    // Need to negate the x and z values since the value is stored for even permutation of the basis vectors,
    // but they are for an odd permutation.  This is an implementation detail.
    /// @brief Return the x-axis coordinate of the point
    inline space::algebra::element_t Point_x(const Point_t &p) {
        return -p.template element<E1_CONF>();
    };

    /// @brief Return the y-axis coordinate of the point
    inline space::algebra::element_t Point_y(const Point_t &p) {
        return p.template element<E2_CONF>();
    };

    /// @brief Return the z-axis coordinate of the point
    inline space::algebra::element_t Point_z(const Point_t &p) {
        return -p.template element<E3_CONF>();
    };

    namespace detail {
        // Support for a Point object.  PGA3 points are sums of trivectors.

        // go through evaluation type checks for the required grade
        // value = 3    : trivector only elements
        // value = -1   : not a point type
//        template<class A>
//        struct grade_evaluation_type {
//            static const int value =
//                    (check_trivector<typename A::clist>::value) ? 3 :
//                    -1;
//        };
//
//        template<class A, int ET = grade_evaluation_type<A>::value>
//        struct Point : public expression<unit<A>> {
//            static_assert(ET != -1, "doesn't evaluate to a PGA3 point");
//        };

        /// @brief Construct a point object

        template<class A>
        struct Point : public gaalet::expression<Point<A>> {
            typedef typename A::clist clist;

            typedef typename A::metric metric;

            typedef typename A::element_t element_t;

            explicit Point(const A &a_)
                    : a(a_) {}

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

            /// @brief Return the multivector representation
            Point_t mv() const {
                return origin() * E0 + x() * E1 + y() * E2 + z() * E3;
            }

            /// @brief Return the normalized multivector representation
            Point_t normalized() const {
                return pga3::normalize(a);
            }

	        /// @brief Return the proper point (where e123 == 1) multivector representation
	        Point_t proper() const {
		        Point_t p = pga3::normalize(a);
		        return (1. / p.template element<E0_CONF>()) * p;
	        }

        protected:
            const A a;
        };

    } // end of namespace detail

    /// @brief PGA3 point
    template <class A> inline
    detail::Point<gaalet::grade<3, A>>
    Point(const gaalet::expression<A>& a) {
        auto g3_a = ::grade<3>(a);
        return detail::Point<decltype(g3_a)>(g3_a);
    }

    /// @brief Point as a multivector object
    template<class A> inline
    gaalet::multivector<typename A::clist, typename A::metric, typename A::element_t>
    eval(const detail::Point<A>& a) {
        return detail::Point<A>(a).eval();
    }
} // end of pga3 namespace

//not allowed, but necessary to work within other scopes ("Argument-dependent name lookup")
namespace std
{
    /// @brief Point object to basic_ostream
    template<class E, class T, class A>
    std::basic_ostream<E, T>& operator<<(std::basic_ostream<E, T>& out, const pga3::detail::Point<A> &p)
    {
        return out << "(x: " << std::fixed << std::setprecision(5) << p.x()
                   << ", y: " << std::fixed << std::setprecision(5) << p.y()
                   << ", z: " << std::fixed << std::setprecision(5) << p.z()
                   << ", O: " << std::fixed << std::setprecision(5) << p.origin()
                   << ")";
    }

} // end of std namespace

#endif // __GAALET_PGA3_POINT_H
