#ifndef __GAALET_PGA3_NORMALIZE_H
#define __GAALET_PGA3_NORMALIZE_H

#pragma once

// Normalize a PGA3 k-vector consisting of a single grade
// A normalized k-vector squares to +/-1, except for some ideal vectors that square to 0

#include "grade_check.h"
#include "magnitude.h"
#include "utility.h"
#include "pga3.h"
#include "pga3_dual.h"

namespace pga3 {
    namespace detail {

        // check for configuration list that describes a PGA3 pseudoscalar only
        template <typename CL>
        struct check_pseudoscalar {
            static const bool value = CL::size == 1 && CL::head == pseudoscalar_conf;
        };

        // check for multivectors of a single grade that can be normalized
        // value=0 : scalar only
        // value=1 : vector only elements
        // value=2 : bivector only elements
        // value=3 : trivector only elements
        // value=4 : pseudoscalar only
        // value=-1 : a mix of grades or an invalid PGA3 grade
        template<class A>
        struct normalize_evaluation_type {
            static const int value =
                    (check_scalar<typename A::clist>::value) ? 0 :
                    (check_vector<typename A::clist>::value) ? 1 :
                    (check_bivector<typename A::clist>::value) ? 2 :
                    (check_trivector<typename A::clist>::value) ? 3 :
                    (check_pseudoscalar<typename A::clist>::value) ? 4 :
                    -1;
        };

        template <class A, int ET = normalize_evaluation_type<A>::value>
        struct unit : public expression<unit<A>> {
            static_assert(ET != -1, "normalization of non-single grade multivectors has not been implemented");
        };

        // Scalar normalization (always return +/- 1)
        template<class A>
        struct unit<A, 0> : public expression<unit<A>> {
            typedef typename A::clist clist;

            typedef typename A::metric metric;

            typedef typename A::element_t element_t;

            unit(const A &x_)
                    : x(x_) {}

            template<conf_t conf>
            element_t element() const {
                return element_t(x.template element<conf>() < 0 ? -1 : 1);
            }
            protected:
                const A x;
        };

        // Vector normalization
        template<class A>
        struct unit<A, 1> : public expression<unit<A>> {
            typedef typename A::clist clist;

            typedef typename A::metric metric;

            typedef typename A::element_t element_t;

            unit(const A &x_)
                    : x(x_)
                    , the_norm(1)
                    ,  first_eval(true) {}

            template<conf_t conf>
            element_t element() const {
                if (first_eval) {
                    the_norm = ::magnitude<A>(x).template element<0>();
                    if (the_norm == 0) {
                        // The vector is ideal, so the norm is just the multiple of e0
                        // Note that this squares to 0
                        the_norm = fabs(x.template element<vector_conf(0)>());
                    }
                    first_eval = false;
                }

                return x.template element<conf>() * (1. / the_norm);
            }
        protected:
            const A x;
            mutable element_t the_norm;
            mutable bool first_eval;
        };

        // Bivector normalization
        template<class A>
        struct unit<A, 2> : public expression<unit<A>> {
            typedef typename A::clist clist;

            typedef typename A::metric metric;

            typedef typename A::element_t element_t;

            unit(const A &x_)
                    : x(x_)
                    , the_norm(1)
                    ,  first_eval(true) {}

            template<conf_t conf>
            element_t element() const {
                if (first_eval) {
                    the_norm = ::magnitude<A>(x).template element<0>();
                    if (the_norm == 0) {
                        // The bivector is ideal, so find the magnitude when joined with the origin (which is
                        // a representative Euclidean point) to compute the ideal norm
                        // Note that e0 is the dual of E0
                        auto x_joined_with_origin = pga3::dual((pga3::dual(x) ^ e0));
                        the_norm = (x_joined_with_origin * (~x_joined_with_origin)).template element<0>();
                    }
                    first_eval = false;
                }

                return x.template element<conf>() * (1. / the_norm);
            }

        protected:
            const A x;
            mutable element_t the_norm;
            mutable bool first_eval;
        };

        // Trivector normalization
        template<class A>
        struct unit<A, 3> : public expression<unit<A>> {
            typedef typename A::clist clist;

            typedef typename A::metric metric;

            typedef typename A::element_t element_t;

            unit(const A &x_)
                    : x(x_)
                    , the_norm(1)
                    ,  first_eval(true) {}

            template<conf_t conf>
            element_t element() const {
                if (first_eval) {
                    the_norm = ::magnitude<A>(x).template element<0>();
                    if (the_norm == 0) {
                        // The trivector (the point) is ideal, so find the magnitude when joined with the origin, E0
                        // (which is a representative Euclidean point) to compute the ideal norm
                        // Note that e0 is the dual of E0
                        auto x_joined_with_origin = pga3::dual((pga3::dual(x) ^ e0));
                        the_norm = (x_joined_with_origin * (~x_joined_with_origin)).template element<0>();
                    }
                    first_eval = false;
                }

                return x.template element<conf>() * (1. / the_norm);
            }

        protected:
            const A x;
            mutable element_t the_norm;
            mutable bool first_eval;
        };

        // Pseudoscalar normalization (always return +/- 1)
        template<class A>
        struct unit<A, 4> : public expression<unit<A>> {
            typedef typename A::clist clist;

            typedef typename A::metric metric;

            typedef typename A::element_t element_t;

            unit(const A &x_)
                    : x(x_) {}

            template<conf_t conf>
            element_t element() const {
                return element_t(x.template element<conf>() < 0 ? -1 : 1);
            }
        protected:
            const A x;
        };

    }  //end namespace detail

    /// @brief Normalize a single grade PGA3 multivector
    /// @ingroup ga_ops
    template<class A>
    inline
    detail::unit<A>
    normalize(const gaalet::expression<A> &a) {
        return detail::unit<A>(a);
    }

} // end of namespace pga3

#endif // __GAALET_PGA3_NORMALIZE_H