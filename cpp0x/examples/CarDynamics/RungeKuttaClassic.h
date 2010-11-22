#ifndef __Runge_Kutta_Classic_h
#define __Runge_Kutta_Classic_h

#include <tuple>

//Wrapper class for CRTP
template <class E>
struct tuple_expression {
   operator const E& () const {
      return *static_cast<const E*>(this);
   }
};

template<typename L, typename R>
struct tuple_addition : public tuple_expression<tuple_addition<L,R>>
{
   tuple_addition(const L& l_, const R& r_)
      :  l(l_),
         r(r_)
   { }

   template<int I>
   auto element() -> decltype(std::get<I>(L()) + std::get<I>(R())) {
   //typename std::tuple_element<I, tuple_t>::type element() {
      return std::get<I>(l) + std::get<I>(r);
   }

   const L& l;
   const R& r;
};

template<typename L>
struct tuple_scalar_product : public tuple_expression<tuple_scalar_product<L>>
{
   tuple_scalar_product(const L& l_, const double& r_)
      :  l(l_),
         r(r_)
   { }

   template<int I>
   auto element() -> decltype(std::get<I>(L())*double()) {
   //typename std::tuple_element<I, L>::type element() {
      return std::get<I>(l)*r;
   }

   const L& l;
   double r;
};

template<class L, class R> inline
tuple_addition<L, R>
operator+(const tuple_expression<L>& l, const tuple_expression<R>& r) {
   return tuple_addition<L, R>(l, r);
}

template<class L, class... A> inline
tuple_addition<L, std::tuple<A...>>
operator+(const tuple_expression<L>& l, const std::tuple<A...>& r) {
   return tuple_addition<L, std::tuple<A...>>(l, r);
}

template<class... A, class R> inline
tuple_addition<std::tuple<A...>, R>
operator+(const std::tuple<A...>& l, const tuple_expression<R>& r) {
   return tuple_addition<std::tuple<A...>, R>(l, r);
}

template<class... L, class... R> inline
tuple_addition<std::tuple<L...>, std::tuple<R...>>
operator+(const std::tuple<L...>& l, const std::tuple<R...>& r) {
   return tuple_addition<std::tuple<L...>, std::tuple<R...>>(l, r);
}

template<class T> inline
tuple_scalar_product<T>
operator*(const tuple_expression<T>& l, const double& r) {
   return tuple_scalar_product<T>(l, r);
}

template<class... A> inline
tuple_scalar_product<std::tuple<A...>>
operator*(const std::tuple<A...>& l, const double& r) {
   return tuple_scalar_product<std::tuple<A...>>(l, r);
}

template<class... A> inline
tuple_scalar_product<std::tuple<A...>>
operator*(const double& l, const std::tuple<A...>& r) {
   return tuple_scalar_product<std::tuple<A...>>(r, l);
}

template<typename E, typename R, int I=std::length<R>::value>
struct tuple_expression_evaluation
{
   static void operator()(const E& e, R& r)
   {
      std::get<I>(r) = e.element<I>();
   }
};

#endif
