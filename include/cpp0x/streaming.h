#ifndef __GAALET_STREAMING_H
#define __GAALET_STREAMING_H

#include <cmath>
#include "expression.h"
#include "multivector.h"

namespace gaalet
{

//expression streaming
template<typename G, typename clist>
struct UnpackElementsToStream
{
   template<class E, class T>
   //static void unpack(std::basic_ostream<E, T>& os, const G& e) {
   static void unpack(std::basic_ostream<E, T>& os, const gaalet::expression<G>& e_) {
      const G& e(e_);
      os << e.template element<clist::head>() << " ";
      UnpackElementsToStream<G, typename clist::tail>::unpack(os, e_);
   }
};

template<typename G>
struct UnpackElementsToStream<G, gaalet::cl_null>
{
   template<class E, class T>
   //static void unpack(std::basic_ostream<E, T>&, const G&) { }
   static void unpack(std::basic_ostream<E, T>&, const gaalet::expression<G>&) { }
};

template<typename clist>
struct UnpackConfigurationListToStream
{
   template<class E, class T>
   static void unpack(std::basic_ostream<E, T>& os) {
      os << clist::head << " ";
      UnpackConfigurationListToStream<typename clist::tail>::unpack(os);
   }
};

template<>
struct UnpackConfigurationListToStream<gaalet::cl_null>
{
   template<class E, class T>
   static void unpack(std::basic_ostream<E, T>&) { }
};

} //end namespace gaalet

//not allowed, but necessary to work within other scopes ("Argument-dependent name lookup")
namespace std
{

template<class E, class T, class G>
std::basic_ostream<E, T>& operator<<(std::basic_ostream<E, T>& os, const gaalet::expression<G>& e)
{
   //const G& e(e_);
   //auto mv = eval(e_);

   os << "[ " << std::dec << std::setprecision(10);
      gaalet::UnpackElementsToStream<G, typename G::clist>::unpack(os, e);
   os << "] { " << std::hex;
      gaalet::UnpackConfigurationListToStream<typename G::clist>::unpack(os);
   os << '}';

   return os;
}

template<class E, class T, class G>
std::basic_ostream<E, T>& operator<<(std::basic_ostream<E, T>&& os, const gaalet::expression<G>& e)
{
   return (os << e);
}

}  // end namespace std

#endif
