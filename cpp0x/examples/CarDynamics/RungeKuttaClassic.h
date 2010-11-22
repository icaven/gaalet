#ifndef __Runge_Kutta_Classic_h
#define __Runge_Kutta_Classic_h

#include <tuple>

template<typename... A>
struct tuple_addition
{
   typedef std::tuple<A...> tuple_t;

   tuple_addition(const tuple_t& l_, const tuple_t& r_)
      :  l(l_),
         r(r_)
   { }

   template<int I>
   //auto element() -> decltype(std::get<I>(l) + std::get<I>(r)) {
   typename std::tuple_element<I, tuple_t>::type element() {
      return std::get<I>(l) + std::get<I>(r);
   }

   const tuple_t& l;
   const tuple_t& r;
};

#endif
