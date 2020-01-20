#pragma once

#include <cmath>

#include "geometric_product.h"
#include "grade.h"
#include "inverse.h"
#include "magnitude.h"

namespace pga3 {
namespace detail {    

//check for general rotor
template<typename CL>
struct check_motor
{
   static const bool value = (BitCount<CL::head>::value==4 || BitCount<CL::head>::value==2 || BitCount<CL::head>::value==0)
                              ? check_motor<typename CL::tail>::value : false;
};
template<>
struct check_motor<cl_null>
{
   static const bool value = true;
};

// check for spinor
template <typename CL>
struct check_spinor {
    static const bool value = (BitCount<CL::head>::value == 2 || BitCount<CL::head>::value == 0) ?
        check_spinor<typename CL::tail>::value :
        false;
};
template <>
struct check_spinor<cl_null> {
    static const bool value = true;
};

// go through inversion evaluation type checks
// value=1 - bivector logarithm
// value=0 - scalar logarithm
// value=2 - spinor logarithm
//value=3 - general rotor logarithm
template <class A>
struct logarithm_evaluation_type {
   static const int value = (check_bivector<typename A::clist>::value) ? 1 :
                            (check_motor<typename A::clist>::value) ? 3 :
                            (check_scalar<typename A::clist>::value) ? 0 :
                            (check_spinor<typename A::clist>::value) ? 2 :
                            -1;
};

template <class A, int ET = logarithm_evaluation_type<A>::value>
struct logarithm : public expression<logarithm<A>> {
    static_assert(ET != -1, "no method for evaluating this type of multivector implemented");
};

//bivector logarithm
template<class A>
struct logarithm<A, 1> : public expression<logarithm<A>>
{
   typedef typename insert_element<0, typename A::clist>::clist clist;

   typedef typename A::metric metric;

   typedef typename A::element_t element_t;

   //dangerous implementation: constructor only called when expression is defined, not when evaluated
   logarithm(const A& a_)
      :  a(a_),
         first_eval(true)
   { }

   //review: don't evaluate on definition workaround: will only work if arguments stay the same (thus attention with variables)
   template<conf_t conf>
   element_t element() const {
      if(first_eval) {
         auto b = eval(::grade<2>(a));
         typedef decltype(b) b_type;
         element_t b_square = 0.0;
         for(int bIt = 0; bIt < b_type::size; ++bIt) {
            b_square += b[bIt]*b[bIt];
         }
         inv_mag_b = 1.0/sqrt(b_square);
         element_t r = a.template element<0>();
         mag_s = sqrt(r*r+b_square);
         first_eval = false;
      }

      return conf==0x00 ? log(mag_s) : a.template element<conf>()*inv_mag_b*acos(0);
   }
   /*template<>
   element_t element<0>() const {
      return ca;
   }*/

protected:
   const A a;
   mutable element_t mag_s;
   mutable element_t inv_mag_b;
   mutable bool first_eval;
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
struct logarithm<A, 3> : public expression<logarithm<A>> {
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
    mutable A B;            // The bivector axis
    mutable A B_polar;      // The axis that is polar to the bivector axis
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

// spinor logarithm
template <class A>
struct logarithm<A, 2> : public expression<logarithm<A>> {
    typedef typename insert_element<0, typename A::clist>::clist clist;

    typedef typename A::metric metric;

    typedef typename A::element_t element_t;

    // dangerous implementation: constructor only called when expression is defined, not when evaluated
    logarithm(const A& a_)
        : a(a_)
        , first_eval(true)
    {
    }

    // review: don't evaluate on definition workaround: will only work if arguments stay the same (thus attention with
    // variables)
    template <conf_t conf>
    element_t element() const
    {
        if(first_eval) {
            // element_t b_square = eval(::magnitude(::grade<2>(a)));
            auto b = eval(::grade<2>(a));
            typedef decltype(b) b_type;
            element_t b_square = 0.0;
            for(unsigned int bIt = 0; bIt < b_type::size; ++bIt) {
                b_square += b[bIt] * b[bIt];
            }
            element_t r = a.template element<0>();
            mag_s = sqrt(r * r + b_square);
            b_acos_r_s = acos(r / mag_s) / sqrt(b_square);
            first_eval = false;
        }

        return conf == 0x00 ? log(mag_s) : a.template element<conf>() * b_acos_r_s;
    }
    /*template<>
    element_t element<0>() const {
       return ca;
    }*/

protected:
    const A a;
    mutable element_t mag_s;
    mutable element_t b_acos_r_s;
    mutable bool first_eval;
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
