#pragma once

#include <cmath>

#include "geometric_product.h"
#include "grade.h"
#include "inverse.h"
#include "magnitude.h"

namespace pga3 { namespace detail {

template <typename CL>
struct check_scalar {
    static const bool value = CL::size == 1 && BitCount<CL::head>::value == 0;
};

//check that all grades in the multivector are even
template<typename CL>
struct check_even_grade
{
    static const bool value = (BitCount<CL::head>::value % 2 == 0) ? check_even_grade<typename CL::tail>::value : false;
};
template<>
struct check_even_grade<cl_null>
{
    static const bool value = true;
};

// go through logarithm evaluation type checks
// value=0 - scalar logarithm
// value=1 - even grade logarithm

template <class A>
struct logarithm_evaluation_type {
    static const int value = (check_scalar<typename A::clist>::value) ? 0 :
                             (check_even_grade<typename A::clist>::value) ? 1 :
                             -1;
};

template <class A, int ET = logarithm_evaluation_type<A>::value>
struct logarithm : public expression<logarithm<A>> {
    static_assert(ET != -1, "no method for evaluating this type of multivector implemented");
};

// motor logarithm
// For PGA, the motor may be expressed as:
// e^((t+uI)B) = s1 + (s2 + p2*I)*B + p1*I
// where I is the pseudoscalar, so I^2 = 0 (for PGA signature 3, 0, 1)
// and t, u are real scalars and B is a bivector that is the axis of the bivector (s2 + p2*I)*B
//
// There are two cases for the axis B: B^2 = 0 or B^2 = -1
// If B^2 == 0
// This means that the exponential is in the form: e^((t+uI)B) = 1 + (t + u*I)*B  (see Gunn2011, Table 7.4)
// and t = s2/s1 and u = p2/s1, unless s1 == 0, in which case u is arbitrary (set it to 0), and set t = 1

//
// However, if B^2 == -1:
//   the screw is of the form: e^((t+uI)B) = cos(t) + (sin (t) + cos (t)uI)B − sin (t)*u*I
// so    s1 = cos(t)
//       s2 = sin(t)
//       p1 = -sin(t) * u
//       p2 = cos(t) * u
//
// The twist, (t+uI)B, is the logarithm of the motor:

template <class A>
struct logarithm<A, 1> : public expression<logarithm<A>> {
    static const conf_t pseudoscalar_conf = Power<2, A::metric::dimension>::value - 1;

   // The output configuration list is composed of all of the bivector elements
   typedef typename insert_element<I_CONF, 
           typename insert_element<J_CONF, 
           typename insert_element<K_CONF, 
           typename insert_element<DUAL_K_CONF, 
           typename insert_element<DUAL_J_CONF, 
           typename insert_element<DUAL_I_CONF, cl_null>::clist>::clist>::clist>::clist>::clist>::clist clist;

    typedef typename A::metric metric;

    typedef typename A::element_t element_t;

    // dangerous implementation: constructor only called when expression is defined, not when evaluated
    logarithm(const A& a_)
        : a(a_)
        , B(pga3::Line_t())
        , B_polar(pga3::Line_t())
        , t(1.0)
        , u(0.0)
        , first_eval(true)
    {
    }

    // review: don't evaluate on definition workaround: will only work if arguments stay the same (thus attention with
    // variables)
    template <conf_t conf>
    element_t element() const
    {
        if(first_eval) {
            element_t s1 = ::grade<0>(a).template element<0>();
            element_t p1 = ::grade<metric::dimension>(a).template element<pseudoscalar_conf>();
            Line_t bivector = ::grade<2>(a);

            auto study_number_and_axis = pga3::bivector_axis(bivector);
            auto z = study_number_and_axis.first;  // study number of the bivector (s2, p2)
            B = study_number_and_axis.second; // axis of the bivector
            B_polar = B * I;
            auto s2 = z.first;
            auto p2 = z.second;
//            std::cout << "s1: " << s1 << " p1: " << p1 << " s2: " << s2 << " p2: " << p2 << std::endl;
//            std::cout << "B: " << B << " B_Polar: " << B_polar << std::endl;

            t = 1.0;
            u = 0.0;

            double B2 = ::part<0>(B * B).template element<0>();
            if(isclose(B2, -1.0)) {
				// Modified from the SIGGRAPH 2019 course notes, p. 44
				if (!isclose(s1, 0.0)) {
					t = atan2(s2, s1);
					u = p2 / s1;
				}
				else if (!isclose(p2, 0.0)) {
					t = atan2(-p1, p2);
					u = -p1 / s2;
				}
			//  else:
			//  both s1 and p2 are zero, the bivector is simple,
			//  u is arbitrary so leave it set to 0, and t to 1

            }
			//else:  # b_squared == 0, therefore  b is an ideal line
				// Translation only - don't modify the default values above
				// t = 1.0
				// u = 0.0
            
            first_eval = false;
        }

        // Result is ((t * pga3::one + u * pga3::I) * B) 
        return  t * B.template element<conf>() + u * B_polar.template element<conf>();
    }

protected:
    const A a;				// The input value
    mutable pga3::Line_t B;            // The bivector axis
    mutable pga3::Line_t B_polar;      // The axis that is polar to the bivector axis
    mutable element_t t;
    mutable element_t u;
    mutable bool first_eval;
};

// scalar logarithm
template <class A>
struct logarithm<A, 0> : public expression<logarithm<A>> {
    typedef typename insert_element<0, cl_null>::clist clist;

    typedef typename A::metric metric;

    typedef typename A::element_t element_t;

    logarithm(const A& a_)
        : a(a_)
    {
    }

    template <conf_t conf>
    element_t element() const
    {
        return (conf == 0) ? log(a.template element<conf>()) : 0.0;
    }

protected:
    const A a;
};

} // end namespace pga3::detail

/// Logarithm of a multivector.
/**
 * Only implemented for scalars, bivectors and spinors.
 */
/// \ingroup ga_ops
template <class A>
inline detail::logarithm<A> log(const gaalet::expression<A>& a)
{
    return detail::logarithm<A>(a);
}

} // end of namespace pga3
