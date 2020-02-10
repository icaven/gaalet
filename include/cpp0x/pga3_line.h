#ifndef __GAALET_PGA3_LINE_H
#define __GAALET_PGA3_LINE_H

#pragma once

#include "grade_check.h"
#include "pga3.h"
#include "pga3_normalize.h"

//
// PGA3 Line support
//

namespace pga3 {

    /// @brief Lines can be defined by Plücker coordinates.  Ensure it is normalized.
    template<typename TX, typename TY, typename TZ, typename DX, typename DY, typename DZ>
    inline
    auto make_line(TX px, TY py, TZ pz, DX dx, DY dy, DZ dz) {
        return pga3::normalize(dx * dual_k + dy * dual_j + dz * dual_i + px * k + py * j + pz * i);
    }

    namespace detail {

        // Support for a Line object.  PGA3 lines are sums of bivectors.


        /// @brief Construct a line object

        template<class A>
        struct Line : public gaalet::expression<Line<A>> {
            typedef typename A::clist clist;

            typedef typename A::metric metric;

            typedef typename A::element_t element_t;

            explicit Line(const A &a_)
                    : a(a_) {}

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

            Line_t eval() const {
                return make_line(k(), j(), i(), dual_k(), dual_j(), dual_i());
            }

            Line_t normalized() const {
                return pga3::normalize(a);
            }

        protected:
            const A a;
        };
    } // end of namespace detail

    /// @brief PGA3 line
    template <class A> inline
    detail::Line<gaalet::grade<2, A>>
    Line(const gaalet::expression<A>& a) {
        auto g2_a = ::grade<2>(a);
        return detail::Line<decltype(g2_a)>(g2_a);
    }

    /// @brief Line as a multivector object
    template<class A> inline
    gaalet::multivector<typename A::clist, typename A::metric, typename A::element_t>
    eval(const detail::Line<A>& a) {
        return detail::Line<A>(a).eval();
    }

} // end of pga3 namespace

//not allowed, but necessary to work within other scopes ("Argument-dependent name lookup")
namespace std
{

    /// @brief Line object to ostream
    template<class A>
    std::ostream &operator<<(std::ostream &out, const pga3::detail::Line<A> &l) {
        return out << "(i: " << l.i()
                   << ", j: " << l.j()
                   << ", k: " << l.k()
                   << ", di: " << l.dual_i()
                   << ", dj: " << l.dual_j()
                   << ", dk: " << l.dual_k()
                   << ")";
    }


} // end of std namespace

#endif // __GAALET_PGA3_LINE_H
