#ifndef __GAALET_INVERSE_H
#define __GAALET_INVERSE_H

#include "part.h"
#include "reverse.h"

namespace gaalet {

}  //end namespace gaalet


template <class A> inline
gaalet::multivector<typename A::clist, A::signature>
operator!(const gaalet::expression<A>& a)
{
   gaalet::element_t div = 1.0/((~a)*a).template element<0x00>();
   return eval(a*div);
}

#endif
