#ifndef __GAALET_PGA3_PLANE_H
#define __GAALET_PGA3_PLANE_H

#pragma once

#include "grade_check.h"
#include "pga3.h"
#include "pga3_normalize.h"

//
// PGA3 Plane support
//

namespace pga3 {
    
    // Four values define a plane
    template<typename T1, typename T2, typename T3, typename T4>
    inline
    auto make_plane(T1 a, T2 b, T3 c, T4 d) {
        return space::algebra::element_t(d) * e0 + space::algebra::element_t(a) * e1 + space::algebra::element_t(b) * e2 + space::algebra::element_t(c) * e3;
    }

    namespace detail {
        // Support for a Plane object.  PGA3 planes are sums of vectors.

        /// @brief Construct a PGA3 plane object

        template<class A>
        struct Plane : public gaalet::expression<Plane<A>> {
            typedef typename A::clist clist;

            typedef typename A::metric metric;

            typedef typename A::element_t element_t;

            explicit Plane(const A &x_)
                    : x(x_) {}

            template<conf_t conf>
            element_t element() const {
                return x.template element<conf>();
            }

            element_t a() const {
                return x.template element<vector_conf(1)>();
            }

            element_t b() const {
                return x.template element<vector_conf(2)>();
            }

            element_t c() const {
                return x.template element<vector_conf(3)>();
            }

            element_t d() const {
                return x.template element<vector_conf(0)>();
            }

            /// @brief Return the multivector representation
            Plane_t mv() const {
                return make_plane(a(), b(), c(), d());
            }

            /// @brief Return the normalized multivector representation
            Plane_t normalized() const {
                return pga3::normalize(x);
            }

        protected:
            const A x;
        };

    } // end of namespace detail

    /// @brief PGA3 plane
    template <class A> inline
    detail::Plane<gaalet::grade<1, A>>
    Plane(const gaalet::expression<A>& a) {
        auto g1_a = ::grade<1>(a);
        return detail::Plane<decltype(g1_a)>(g1_a);
    }

    /// @brief Plane as a multivector object
    template<class A> inline
    gaalet::multivector<typename A::clist, typename A::metric, typename A::element_t>
    eval(const detail::Plane<A>& a) {
        return detail::Plane<A>(a).eval();
    }

} // end of namespace pga3

//not allowed, but necessary to work within other scopes ("Argument-dependent name lookup")
namespace std
{
    /// @brief Plane object to ostream
    template<class A>
    std::ostream &operator<<(std::ostream &out, const pga3::detail::Plane<A> &p) {
        return out << "(a: " << p.a()
                   << ", b: " << p.b()
                   << ", c: " << p.c()
                   << ", d: " << p.d()
                   << ")";
    }

} // end of std namespace

#endif // __GAALET_PGA3_PLANE_H
