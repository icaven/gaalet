#ifndef __GAALET_PGA3_NORM_H
#define __GAALET_PGA3_NORM_H

#pragma once

// Compute the norm of a PGA3 multivector

#include "grade_check.h"
#include "utility.h"
#include "pga3.h"
#include "pga3_dual.h"
#include "pga3_utility.h"

namespace pga3 {
    namespace detail {

        // check for multivectors of that can be normalized
        // value=0 : scalar only
        // value=4 : pseudoscalar only
        // value=-1 : vector, bivector, trivector or a mix of grades
        template<class A>
        struct norm_evaluation_type {
            static const int value =
                    (check_scalar<typename A::clist>::value) ? 0 :
                    (check_quadvector<typename A::clist>::value) ? 4 :
                    -1;
        };

        template <class A, int ET = norm_evaluation_type<A>::value>
        struct norm2 : public expression<norm2<A>> {
            typedef typename A::clist clist;

            typedef typename A::metric metric;

            typedef typename A::element_t element_t;

            explicit norm2 (const A& x_)
                : m_x(x_)
                , m_first_eval(true)
                {};
            
            template <conf_t conf>
            element_t element() {
                if (m_first_eval) {
                    m_norm2 = ::grade<0>(m_x*(~m_x)).template element<0>();
                    if (pga3::isclose(m_norm2, 0)) {
                        // The element is ideal, so find the magnitude when joined with the origin, E0
                        // (which is a representative Euclidean point) to compute the ideal norm
                        // Note that e0 is the dual of E0
                        auto x_joined_with_origin = pga3::dual((pga3::dual(m_x) ^ e0));
                        m_norm2 = ::grade<0>(x_joined_with_origin * (~x_joined_with_origin)).template element<0>();
                    }
                    m_first_eval = false;
                }
                return (conf==0x00) ? m_norm2 : element_t(0.0);
            }
            
        protected:
            const A m_x;
            mutable bool m_first_eval;
            mutable element_t m_norm2;

        };

        // Squared value of scalar
        template<class A>
        struct norm2<A, 0> : public expression<norm2<A>> {
            typedef typename A::clist clist;

            typedef typename A::metric metric;

            typedef typename A::element_t element_t;

            explicit norm2(const A &x_)
                    : x(x_) {}

            template<conf_t conf>
            element_t element() const {
                return x.template element<conf>() * x.template element<conf>();
            }
            protected:
                const A x;
        };
        
        // Pseudoscalar squared norm
        template<class A>
        struct norm2<A, 4> : public expression<norm2<A>> {
            typedef typename A::clist clist;

            typedef typename A::metric metric;

            typedef typename A::element_t element_t;

            explicit norm2(const A &x_)
                    : x(x_) {}

            template<conf_t conf>
            element_t element() const {
                return x.template element<conf>() * x.template element<conf>();
            }
        protected:
            const A x;
        };


    }  //end namespace detail

    /// @brief Compute the squared norm of a PGA3 multivector
    /// @ingroup ga_ops
    template<class A>
    inline
    typename A::element_t
    norm2(const gaalet::expression<A> &a) {
        return std::max(detail::norm2<A>(a).template element<0>(),
                        detail::norm2<A>(a).template element<pseudoscalar_conf>());
    }

    /// @brief Compute the norm of a PGA3 multivector
    /// @ingroup ga_ops
    template<class A>
    inline
    typename A::element_t
    norm(const gaalet::expression<A> &a) {
        typename A::element_t squared_magnitude = std::max(detail::norm2<A>(a).template element<0>(),
                                                           detail::norm2<A>(a).template element<pseudoscalar_conf>());
        return sqrt(squared_magnitude);
    }
} // end of namespace pga3

#endif //__GAALET_PGA3_NORM_H