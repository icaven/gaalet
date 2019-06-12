#pragma once

#include <cmath>
#include <utility>

#include "geometric_product.h"
#include "grade.h"

#include "pga3.h"

//
// Specialized version of the exponential function for the PGA algebra
// TODO: modify the library exponential function to support this algebra, since this version
// is only a slight modification of the library function

namespace pga3
{
namespace detail
{

    // check for bivector
    template <typename CL>
    struct check_bivector {
        static const bool value = (BitCount<CL::head>::value == 2) ? check_bivector<typename CL::tail>::value : false;
    };
    template <>
    struct check_bivector<cl_null> {
        static const bool value = true;
    };

    // check for scalar
    template <typename CL>
    struct check_scalar {
        static const bool value = (CL::size == 1 && BitCount<CL::head>::value == 0) ? true : false;
    };

    // go through inversion evaluation type checks
    // value=1 - bivector exponential
    // value=2 - scalar exponential
    template <class A>
    struct exponential_evaluation_type {
        static const int value =
            (check_bivector<typename A::clist>::value) ? 1 : (check_scalar<typename A::clist>::value) ? 0 : -1;
    };

    template <class A, int ET = exponential_evaluation_type<A>::value>
    struct exponential : public expression<exponential<A>> {
        static_assert(ET != -1, "no method for evaluating this type of multivector implemented");
    };

    // bivector exponential
    template <class A>
    struct exponential<A, 1> : public expression<exponential<A>> {
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
                auto a = z.first;
                b = z.second;

//                std::cout << "<sn>:" << a << " " << b << " </sn>" << std::endl;
//                std::cout << "<axis>:" << B << " </axis>" << std::endl;
//                std::cout << "<B_polar>:" << B_polar << " </B_polar>" << std::endl;
//                std::cout << "<a*a>:" << (a * a) << " </a*a>" << std::endl;

                if(a > 0.0) {
                    ca = cos(a);
                    sa = sin(a);
                    ca_b = ca * b;
                    sign_ca = signbit(ca) ? -1.0 : 1.0;
//                    std::cout << "ca:" << ca << " sa: " << sa << " ca_b: " << ca_b << std::endl;
                } else {
                    ca = 1.0;
                    sa = 1.0;
                    ca_b = 0.0;
                    sign_ca = 1.0;
                }
                first_eval = false;
            }

            if(conf == pseudoscalar_conf) {
                return -sign_ca * sa * b;
            } else if(conf == 0) {
                return sign_ca * ca;
            } else {
                return sign_ca * (B.template element<conf>() * sa + B_polar.template element<conf>() * ca_b);
            }
        }

    protected:
        const A a;
        mutable Line_t B;       // The bivector axis
        mutable Line_t B_polar; // The axis that is polar to the bivector axis
        mutable element_t b;
        mutable element_t ca;
        mutable element_t sa;
        mutable element_t ca_b;
        mutable element_t sign_ca;
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
#undef exp

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
