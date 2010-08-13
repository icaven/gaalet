#ifndef __GAALET_EXPRESSION_H
#define __GAALET_EXPRESSION_H

#include "configuration_list.h"


#include <iostream>
#include <iomanip>


namespace gaalet
{

//Wrapper class for CRTP
template <class E>
struct expression {
   operator const E& () const {
      return *static_cast<const E*>(this);
   }
};

} //end namespace gaalet

//expression streaming
template<typename E, typename clist>
struct UnpackElementsToStream
{
   static void unpack(std::ostream& os, const E& e) {
      os << e.template element<clist::head>() << ' ';
      UnpackElementsToStream<E, typename clist::tail>::unpack(os, e);
   }
};

template<typename E>
struct UnpackElementsToStream<E, gaalet::cl_null>
{
   static void unpack(std::ostream&, const E&) { }
};

template<typename clist>
struct UnpackConfigurationListToStream
{
   static void unpack(std::ostream& os) {
      os << clist::head << ' ';
      UnpackConfigurationListToStream<typename clist::tail>::unpack(os);
   }
};

template<>
struct UnpackConfigurationListToStream<gaalet::cl_null>
{
   static void unpack(std::ostream&) { }
};

template<class E>
std::ostream& operator<<(std::ostream& os, const gaalet::expression<E>& e_)
{
   const E& e(e_);

   os << "[ ";
      UnpackElementsToStream<E, typename E::clist>::unpack(os, e);
   os << "] { " << std::hex;
      UnpackConfigurationListToStream<typename E::clist>::unpack(os);
   os << '}';

   return os;
}


#endif
