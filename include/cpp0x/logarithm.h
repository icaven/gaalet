#ifndef __GAALET_LOGARITHM_H
#define __GAALET_LOGARITHM_H

#include <cmath>

#include "grade.h"
#include "geometric_product.h"

namespace gaalet {

//check for bivector
template<typename CL>
struct check_bivector
{
   static const bool value = (BitCount<CL::head>::value==2) ? check_bivector<typename CL::tail>::value : false;
};
template<>
struct check_bivector<cl_null>
{
   static const bool value = true;
};

//check for scalar
template<typename CL>
struct check_scalar
{
   static const bool value = (CL::size==1 && BitCount<CL::head>::value==0) ? true : false;
};

//check for spinor
template<typename CL>
struct check_spinor
{
   static const bool value = check_scalar<CL>::value && check_bivector<CL>::value;
}

//go through inversion evaluation type checks
//value=1 - bivector logarithm
//value=0 - scalar logarithm
//value=2 - spinor logarithm
template<class A>
struct logarithm_evaluation_type
{
   static const int value = (check_bivector<typename A::clist>::value) ? 1 :
                            (check_scalar<typename A::clist>::value) ? 0 :
                            (check_spinor<typename A::clist>::value) ? 2 :
                            -1;
};

template<class A, int ET = logarithm_evaluation_type<A>::value>
struct logarithm : public expression<logarithm<A>>
{
   static_assert(ET!=-1, "no method for evaluating this type of multivector implemented");
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
         inv_mag = eval(!magnitude(a));
         first_eval = false;
      }

      return a.element<conf>()*inv_mag*acos(0);
   }
   /*template<>
   element_t element<0>() const {
      return ca;
   }*/

protected:
   const A& a;
   mutable element_t inv_mag;
   mutable bool first_eval;
};

//scalar logarithm
template<class A>
struct logarithm<A, 0> : public expression<logarithm<A>>
{
   typedef typename insert_element<0, cl_null>::clist clist;

   typedef typename A::metric metric;

   typedef typename A::element_t element_t;

   logarithm(const A& a_)
      :  a(a_)
   { }

   template<conf_t conf>
   element_t element() const {
      return (conf==0) ? exp(a.element<conf>()) : 0.0;
   }

protected:
   const A& a;
};


}  //end namespace gaalet

template <class A> inline
gaalet::logarithm<A>
exp(const gaalet::expression<A>& a) {
   return gaalet::logarithm<A>(a);
}


#endif
