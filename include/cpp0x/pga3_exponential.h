#ifndef __GAALET_PGA3_EXPONENTIAL_H
#define __GAALET_PGA3_EXPONENTIAL_H

#pragma once

#include <cmath>
#include <utility>

#include "geometric_product.h"
#include "grade_check.h"
#include "grade.h"

#include "pga3.h"
#include "pga3_utility.h"

//
// Specialized version of the exponential function for the PGA algebra; see Section 7.8 in [Gunn2011]
//

namespace pga3
{
namespace detail
{
    // go through exponent evaluation type checks
    // value=0 - scalar exponential
    // value=2 - bivector exponential
    template <class A>
    struct exponential_evaluation_type {
        static const int value =
                (check_scalar<typename A::clist>::value) ? 0 :
                (check_bivector<typename A::clist>::value) ? 2  : -1;
    };

    template <class A, int ET = exponential_evaluation_type<A>::value>
    struct exponential : public expression<exponential<A>> {
        static_assert(ET != -1, "no method for evaluating this type of multivector implemented");
    };

    // PGA3 exponential:bivectors
    template <class A>
    struct exponential<A, 2> : public expression<exponential<A>> {
        static const conf_t pseudoscalar_conf = Power<2, A::metric::dimension>::value - 1;

        // The output configuration list is composed of all of the bivector elements and the scalar and pseudoscalar
        typedef typename insert_element<0,
                typename insert_element<I_CONF, 
                typename insert_element<J_CONF, 
                typename insert_element<K_CONF, 
                typename insert_element<DUAL_K_CONF, 
                typename insert_element<DUAL_J_CONF, 
                typename insert_element<DUAL_I_CONF, 
                typename insert_element<pseudoscalar_conf, 
                cl_null>::clist>::clist>::clist>::clist>::clist>::clist>::clist>::clist clist;

        typedef typename A::metric metric;

        typedef typename A::element_t element_t;

        // dangerous implementation: constructor only called when expression is defined, not when evaluated
        exponential(const A& a_)
                : a(a_)
                , B(Line_t())
                , B_polar(Line_t())
                , B2(0.0)
                , f1(1.0)
                , f2(0.0)
                , f3(1.0)
                , f4(1.0)
                , first_eval(true)
        {
        }

        // review: don't evaluate on definition workaround: will only work if arguments stay the same (thus attention
        // with variables)
        template <conf_t conf>
        element_t element() const
        {
            if(first_eval) {
                // Decompose the input bivector into it's Study number and simple bivector axis
                std::pair<std::pair<double, double>, Line_t> study_number_and_axis = bivector_axis(a);
                auto z = study_number_and_axis.first;
                B = study_number_and_axis.second; // axis of the bivector
                B_polar = B * I;
				
				// Specialize this exponential for PGA3D where I*I == 0
				B2 = ::part<0>(B * B).template element<0>();
				if(isclose(B2, -1.0)) {
					f1 = cos(z.first);
					f2 = sin(z.first);
					f3 = 1.0;
					f4 = z.second;
				}
				else { 
					// B*B == 0
					f1 = 1.0;
					f2 = isclose(z.first, 0.) && isclose(z.second, 0.) ? 1.0 : z.first;
					f3 = 1.0;
					f4 = isclose(z.first, 0.) && isclose(z.second, 0.) ? 0.0 : z.second;
				}
			
                first_eval = false;
            }

			return ((f1 * f3) * one +
					(B2 * f2 * f4) * I +
					(f2 * f3) * B + 
					(f1 * f4) * B_polar).template element<conf>();
        }

    protected:
        const A a;
        mutable Line_t B;       // The bivector axis
        mutable Line_t B_polar; // The axis that is polar to the bivector axis
        mutable element_t B2;	// The square of the simple bivector (either is -1 or 0)
		mutable element_t f1;
		mutable element_t f2;
		mutable element_t f3;
		mutable element_t f4;
        mutable bool first_eval;
    };

    // scalar exponential
    template <class A>
    struct exponential<A, 0> : public expression<exponential<A>> {
        typedef typename insert_element<0, cl_null>::clist clist;

        typedef typename A::metric metric;

        typedef typename A::element_t element_t;

        exponential(const A& a_)
            : a(a_)
        {
        }

        template <conf_t conf>
        element_t element() const
        {
            return (conf == 0) ? exp(a.template element<conf>()) : 0.0;
        }

    protected:
        const A a;
    };

} // end namespace detail

// XXX Fix exp() in gaalet, so that this next line isn't needed
//#undef exp

/// Exponential of a multivector.
/**
 * @brief Specialized exponential function for PGA3
 * @param a     The gaalet expression
 * @return      The exponential of a
 *
 * Only implemented for scalars and bivectors.
 */
/// \ingroup ga_ops
template <class A>
inline detail::exponential<A> exp(const gaalet::expression<A>& a)
{
    return detail::exponential<A>(a);
}

} // end namespace pga3

#endif // __GAALET_PGA3_EXPONENTIAL_H
