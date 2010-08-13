#ifndef __GAALET_FAST_MULTIVECTOR_H
#define __GAALET_FAST_MULTIVECTOR_H

#include "configuration_list.h"
#include "signature_list.h"
#include "expression.h"

#include <algorithm>


namespace gaalet
{

//multivector coefficients type
typedef double element_t;

//multivector struct
//template<typename CL, typename SL=sl::sl_null>
//struct multivector : public expression<multivector<CL, SL>>
template<unsigned int ID, typename CL, conf_t S = 0x00>
struct fast_multivector : public expression<fast_multivector<ID, CL, S>>
{
   //fast_multivector id
   static const unsigned int id = ID;

   typedef CL clist;
   static const conf_t size = clist::size;
   
   //typedef SL slist;
   static const conf_t signature = S;

   //initialization
   fast_multivector()
   {
      if(!data) {
         data = new element_t[size];
      }
      std::fill(data, data+size, 0.0);
   }

   fast_multivector(std::initializer_list<element_t> s)
   {
      if(!data) {
         data = new element_t[size];
      }
      element_t* last = std::copy(s.begin(), (s.size()<=size) ? s.end() : (s.begin()+size), data);
      std::fill(last, data+size, 0.0);
   }

   ~fast_multivector()
   {
      delete[] data;
   }

   //return element by index, index known at runtime
   const element_t& operator[](const conf_t& index) const {
      return data[index];
   }
   element_t& operator[](const conf_t& index) {
      return data[index];
   }

   //return element by index, index known at compile time
   template<conf_t index>
   const element_t& get() const {
      return data[index];
   }

   //static return function for fast expression templates
   template<conf_t conf>
   static element_t fast_element() {
      static const conf_t index = search_element<conf, clist>::index;
      //static_assert(index<size, "element<conf_t>(): no such element in configuration list");
      return (index<size) ? data[index] : 0.0;
   }

   //evaluation
   template<typename E, conf_t index = 0>
   struct ElementEvaluation
   {                          //v no reference to pointer *& with gcc4.5 possible... What's going on?
      static void eval(element_t* const data) {
         data[index] = E::template fast_element<get_element<index, clist>::value>();
         ElementEvaluation<E, index+1>::eval(data);
      }
   };
   template<typename E>
   struct ElementEvaluation<E, size-1>
   {
      static void eval(element_t* const data) {
         data[size-1] = E::template fast_element<get_element<size-1, clist>::value>();
      }
   };

   //   constructor evaluation
   template<class E>
   fast_multivector(const expression<E>& e_) {
      if(!data) {
         data = new element_t[size];
      }
      const E& e(e_);
      ElementEvaluation<E>::eval(data, e);
   }

   //   assignment evaluation
   template<class E>
   void operator=(const expression<E>& e_) {
      const E& e(e_);
      ElementEvaluation<E>::eval(data);
   }

protected:
   static element_t* data;
};

template<unsigned int ID, typename CL, conf_t S> element_t* fast_multivector<ID, CL, S>::data = NULL;


//fast_multivector configuration elements unpacking
template<conf_t... elements>
struct fast_mv;

//no cpp0x template aliases supported by gcc yet
/*template<conf_t head, conf_t... tail>
using mv = fast_multivector<typename insert_element<head, typename mv<tail...>::clist>::clist>;

template<>
using mv = fast_multivector<cl_null>;*/

template<conf_t head, conf_t... tail>
struct fast_mv<head, tail...>
{
   typedef fast_multivector<head, typename insert_element<head, typename fast_mv<tail...>::type::clist>::clist> type;
};
template<conf_t head>
struct fast_mv<head>
{
   typedef fast_multivector<head, cl_null> type;
};


} //end namespace gaalet

#endif
