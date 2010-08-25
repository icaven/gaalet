#ifndef __GAALET_EXPRESSION_H
#define __GAALET_EXPRESSION_H

#include "configuration_list.h"


#include <iostream>
#include <iomanip>


namespace gaalet
{

//multivector coefficients type, set hard
//---float type for cuda
#ifdef __CUDACC__
   typedef float element_t;
#else
   typedef double element_t;
#endif

//Wrapper class for CRTP
template <class E>
struct expression {
   operator const E& () const {
      return *static_cast<const E*>(this);
   }
};

} //end namespace gaalet


#endif
