#ifndef __GAALET_UNIT_H
#define __GAALET_UNIT_H

// Normalize a blade

#include "magnitude.h"
#include "utility.h"

namespace gaalet {

// A blade is a outer product of vectors, so if the configuration list is of length 1 
// and if the bitcount of entry isn't 0 (which would indicate a scalar), then the entity is a blade.
template<typename CL>
struct check_blade
{
   static const bool value = CL::size == 1 && BitCount<CL::head>::value > 0;
};

//go through blade evaluation type checks
// value=1:  blade
// value=-1: a general multivector or scalar
template<class A>
struct blade_evaluation_type
{
   static const int value = (check_blade<typename A::clist>::value) ? 1 : -1;
};

//template<class A, int ET=blade_evaluation_type<A>::value>
//struct unit : public expression<unit<A>>
//{
//   static_assert(ET!=-1, "no method for normalizing this type of multivector implemented");
//};

//template<class A>
//struct unit<A, 1> : public expression<unit<A>>
template<class A>
struct unit : public expression<unit<A>>
{    
    typedef typename A::clist clist;

    typedef typename A::metric metric;

    typedef typename A::element_t element_t;

   unit(const A& x_)
      :  x(x_), 
        first_eval(true)
   { }

   template<conf_t conf>
   element_t element() const {
      //review: don't evaluate on definition workaround: will only work if arguments stay the same (thus attention with variables)
        if (conf == (Power<2, A::metric::dimension>::value - 1)) {
            // The normalized pseudoscalar is always +/- 1
            element_t pseudoscalar_multiple = x.template element<conf>();
            int sign = A::metric::q % 2 == 1 ? -1 : 1;
            return pseudoscalar_multiple > 0 ? sign : (pseudoscalar_multiple < 0 ? -sign : 0);
        }

        if (first_eval) {
            the_norm = magnitude<A>(x). template element<0x00>();
            first_eval=false;
        }
        
        if (the_norm != 0) {
            return x.template element<conf>() * (1./ the_norm);
        }
        else {
            return element_t(0.0);                
        }
   }

protected:
   const A x;
   mutable element_t the_norm;
   mutable bool first_eval;
};

}  //end namespace gaalet

/// Unit blade
/**
 * Only implemented for blades. Compile error for general multivectors.
 */
/// \ingroup ga_ops
template <class A> inline
gaalet::unit<A>
normalize(const gaalet::expression<A>& a) {
   return gaalet::unit<A>(a);
}

#endif // __GAALET_UNIT_H