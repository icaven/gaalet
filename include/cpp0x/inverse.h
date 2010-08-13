#ifndef __GAALET_INVERSE_H
#define __GAALET_INVERSE_H

#include "part.h"
#include "reverse.h"

namespace gaalet {

}  //end namespace gaalet

/*template <class A> inline
auto operator!(const gaalet::expression<A>& a) -> decltype((~a)*(1.0/(a*(~a)).template element<0x00>()))
{
   return (~a)*(1.0/(a*(~a)).template element<0x00>());
}*/

template <class A> inline
auto operator!(const gaalet::expression<A>& a) -> decltype(eval(a*gaalet::element_t()))
{
   gaalet::element_t div = 1.0/((~a)*a).template element<0x00>();
   return eval(a*div);
}

#endif
