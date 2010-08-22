#ifndef __GAALET_MAGNITUDE_H
#define __GAALET_MAGNITUDE_H

#include "grade.h"
#include "reverse.h"

template <class A> inline
auto magnitude(const gaalet::expression<A>& a) -> decltype(grade<0>((~a)*a))
{
   return grade<0>((~a)*a);
}

#endif
