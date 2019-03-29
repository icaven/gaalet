#ifndef __GAALET_MAGNITUDE_H
#define __GAALET_MAGNITUDE_H

#include "geometric_product.h"
#include "reverse.h"
#include "grade.h"
#include "symbex.h"  //  This is included only to resolve sqrt() correctly for this type - there should be a better way

namespace gaalet {

template<class A>
struct magnitude : public expression<magnitude<A>>
{
   typedef configuration_list<0x00, cl_null> clist;

   typedef typename A::metric metric;
   
   typedef typename A::element_t element_t;

   magnitude(const A& a_)
      :  a(a_)
   { }

   template<conf_t conf>
   element_t element() const {
      return (conf==0x00) ? ::sqrt(::grade<0>(eval(a*(~a))).template element<0>()) : element_t(0.0);
   }

protected:
   const A a;
};

// The squared magnitude
template<class A>
struct magnitude2 : public expression<magnitude2<A>>
{
   typedef configuration_list<0x00, cl_null> clist;

   typedef typename A::metric metric;
   
   typedef typename A::element_t element_t;

   magnitude2(const A& a_)
      :  a(a_)
   { }

   template<conf_t conf>
   element_t element() const {
      return (conf==0x00) ? ::grade<0>(eval(a*(~a))).template element<0>() : element_t(0.0);
   }

protected:
   const A a;
};

}  //end namespace gaalet

/// \brief Magnitude of a multivector.
/**
 * Following Hestenes' definition. Undefined for degenerate algebra.
 */
/// \ingroup ga_ops
template <class A> inline
gaalet::magnitude<A>
magnitude(const gaalet::expression<A>& a)
{
   return gaalet::magnitude<A>(a);
}

/// \brief Squared magnitude of a multivector.
/**
 * Following Hestenes' definition. Undefined for degenerate algebra.
 */
/// \ingroup ga_ops
template <class A> inline
gaalet::magnitude2<A>
magnitude2(const gaalet::expression<A>& a)
{
   return gaalet::magnitude2<A>(a);
}


#endif
