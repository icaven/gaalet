#ifndef __GAALET_MULTIVECTOR_H
#define __GAALET_MULTIVECTOR_H

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
template<typename CL, conf_t S = 0x00>
struct multivector : public expression<multivector<CL, S> >
{
   typedef CL clist;
   static const conf_t size = clist::size;
   
   //typedef SL slist;
   static const conf_t signature = S;

   //initialization
   multivector()
   {
      std::fill(data, data+size, 0.0);
   }

   multivector(const element_t& e_)
   {
      std::fill(data+1, data+size, 0.0);
      data[0] = e_;
   }


   /*multivector(std::initializer_list<element_t> s)
   {
      element_t* last = std::copy(s.begin(), (s.size()<=size) ? s.end() : (s.begin()+size), data);
      std::fill(last, data+size, 0.0);
   }*/

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

   //return element by configuration (basis vector), configuration known at compile time
   template<conf_t conf>
   ////reference return (const element_t& element() const) not applicable because of possible return of 0.0;
   element_t element() const {
      static const conf_t index = search_element<conf, clist>::index;
      //static_assert(index<size, "element<conf_t>(): no such element in configuration list");
      return (index<size) ? data[index] : 0.0;
   }

   //evaluation
   template<typename E, conf_t index = 0>
   struct ElementEvaluation
   {                          //v no reference to pointer *& with gcc4.5 possible... What's going on?
      static void eval(element_t* const data, const E& e) {
         data[index] = e.element<get_element<index, clist>::value>();
         ElementEvaluation<E, index+1>::eval(data, e);
      }
   };
   template<typename E>
   struct ElementEvaluation<E, size-1>
   {
      static void eval(element_t* const data, const E& e) {
         data[size-1] = e.element<get_element<size-1, clist>::value>();
      }
   };

   //   constructor evaluation
   template<class E>
   multivector(const expression<E>& e_) {
      const E& e(e_);
      ElementEvaluation<E>::eval(data, e);
   }

   //   assignment evaluation
   template<class E>
   void operator=(const expression<E>& e_) {
      const E& e(e_);
      ElementEvaluation<E>::eval(data, e);
   }

protected:
   element_t data[size];
};

//specialization for scalar multivector type
template<conf_t S>
struct multivector<configuration_list<0x00, cl_null>, S> : public expression<multivector<configuration_list<0x00, cl_null>, S> >
{
   typedef configuration_list<0x00, cl_null> clist;
   static const conf_t size = clist::size;
   
   //typedef SL slist;
   static const conf_t signature = S;

   //initialization
   multivector()
      :  value(0.0)
   { }

   multivector(const element_t& setValue)
      :  value(setValue)
   { }

   /*multivector(std::initializer_list<element_t> s)
      :  value(*s.begin())
   { }*/

   operator element_t()
   {
      return value;
   }

   //return element by index, index known at runtime
   const element_t& operator[](const conf_t&) const {
      return value;
   }
   element_t& operator[](const conf_t& index) {
      return value;
   }

   //return element by index, index known at compile time
   template<conf_t index>
   const element_t& get() const {
      return value;
   }

   //return element by configuration (basis vector), configuration known at compile time
   template<conf_t conf>
   ////reference return (const element_t& element() const) not applicable because of possible return of 0.0;
   element_t element() const {
      //static const conf_t index = search_element<conf, clist>::index;
      //static_assert(index<size, "element<conf_t>(): no such element in configuration list");
      return (conf==0x00) ? value : 0.0;
   }

   //   constructor evaluation
   template<class E>
   multivector(const expression<E>& e_) {
      const E& e(e_);
      value = e.element<0x00>();
   }

   //   assignment evaluation
   template<class E>
   void operator=(const expression<E>& e_) {
      const E& e(e_);
      value = e.element<0x00>();
   }

protected:
   element_t value;
};


//multivector configuration elements unpacking
//template<conf_t... elements>
//struct mv;

//no cpp0x template aliases supported by gcc yet
/*template<conf_t head, conf_t... tail>
using mv = multivector<typename insert_element<head, typename mv<tail...>::clist>::clist>;

template<>
using mv = multivector<cl_null>;*/

/*template<conf_t head, conf_t... tail>
struct mv<head, tail...>
{
   typedef multivector<typename insert_element<head, typename mv<tail...>::type::clist>::clist> type;
};
template<>
struct mv<>
{
   typedef multivector<cl_null> type;
};*/

template<conf_t S>
struct metric
{
   static const conf_t signature = S;

   template<conf_t e1=-1, conf_t e2=-1, conf_t e3=-1, conf_t e4=-1, conf_t e5=-1, conf_t e6=-1, conf_t e7=-1, conf_t e8=-1, conf_t e9=-1, conf_t e10=-1, conf_t e11=-1, conf_t e12=-1>
   struct mv
   {
      typedef multivector<typename insert_element<e1,
                          typename insert_element<e2,
                          typename insert_element<e3,
                          typename insert_element<e4,
                          typename insert_element<e5,
                          typename insert_element<e6,
                          typename insert_element<e7,
                          typename insert_element<e8,
                          typename insert_element<e9,
                          typename insert_element<e10,
                          typename insert_element<e11,
                          typename insert_element<e12,
                          cl_null>::clist>::clist>::clist>::clist>::clist>::clist>::clist>::clist>::clist>::clist>::clist>::clist, signature> type;
   };
   template<conf_t e1, conf_t e2, conf_t e3, conf_t e4, conf_t e5, conf_t e6, conf_t e7, conf_t e8, conf_t e9, conf_t e10, conf_t e11>
   struct mv<e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, -1>
   {
      typedef multivector<typename insert_element<e1,
                          typename insert_element<e2,
                          typename insert_element<e3,
                          typename insert_element<e4,
                          typename insert_element<e5,
                          typename insert_element<e6,
                          typename insert_element<e7,
                          typename insert_element<e8,
                          typename insert_element<e9,
                          typename insert_element<e10,
                          typename insert_element<e11,
                          cl_null>::clist>::clist>::clist>::clist>::clist>::clist>::clist>::clist>::clist>::clist>::clist, signature> type;
   };
   template<conf_t e1, conf_t e2, conf_t e3, conf_t e4, conf_t e5, conf_t e6, conf_t e7, conf_t e8, conf_t e9, conf_t e10>
   struct mv<e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, -1, -1>
   {
      typedef multivector<typename insert_element<e1,
                          typename insert_element<e2,
                          typename insert_element<e3,
                          typename insert_element<e4,
                          typename insert_element<e5,
                          typename insert_element<e6,
                          typename insert_element<e7,
                          typename insert_element<e8,
                          typename insert_element<e9,
                          typename insert_element<e10,
                          cl_null>::clist>::clist>::clist>::clist>::clist>::clist>::clist>::clist>::clist>::clist, signature> type;
   };
   template<conf_t e1, conf_t e2, conf_t e3, conf_t e4, conf_t e5, conf_t e6, conf_t e7, conf_t e8, conf_t e9>
   struct mv<e1, e2, e3, e4, e5, e6, e7, e8, e9, -1, -1, -1>
   {
      typedef multivector<typename insert_element<e1,
                          typename insert_element<e2,
                          typename insert_element<e3,
                          typename insert_element<e4,
                          typename insert_element<e5,
                          typename insert_element<e6,
                          typename insert_element<e7,
                          typename insert_element<e8,
                          typename insert_element<e9,
                          cl_null>::clist>::clist>::clist>::clist>::clist>::clist>::clist>::clist>::clist, signature> type;
   };
   template<conf_t e1, conf_t e2, conf_t e3, conf_t e4, conf_t e5, conf_t e6, conf_t e7, conf_t e8>
   struct mv<e1, e2, e3, e4, e5, e6, e7, e8, -1, -1, -1, -1>
   {
      typedef multivector<typename insert_element<e1,
                          typename insert_element<e2,
                          typename insert_element<e3,
                          typename insert_element<e4,
                          typename insert_element<e5,
                          typename insert_element<e6,
                          typename insert_element<e7,
                          typename insert_element<e8,
                          cl_null>::clist>::clist>::clist>::clist>::clist>::clist>::clist>::clist, signature> type;
   };
   template<conf_t e1, conf_t e2, conf_t e3, conf_t e4, conf_t e5, conf_t e6, conf_t e7>
   struct mv<e1, e2, e3, e4, e5, e6, e7, -1, -1, -1, -1, -1>
   {
      typedef multivector<typename insert_element<e1,
                          typename insert_element<e2,
                          typename insert_element<e3,
                          typename insert_element<e4,
                          typename insert_element<e5,
                          typename insert_element<e6,
                          typename insert_element<e7,
                          cl_null>::clist>::clist>::clist>::clist>::clist>::clist>::clist, signature> type;
   };
   template<conf_t e1, conf_t e2, conf_t e3, conf_t e4, conf_t e5, conf_t e6>
   struct mv<e1, e2, e3, e4, e5, e6, -1, -1, -1, -1, -1, -1>
   {
      typedef multivector<typename insert_element<e1,
                          typename insert_element<e2,
                          typename insert_element<e3,
                          typename insert_element<e4,
                          typename insert_element<e5,
                          typename insert_element<e6,
                          cl_null>::clist>::clist>::clist>::clist>::clist>::clist, signature> type;
   };
   template<conf_t e1, conf_t e2, conf_t e3, conf_t e4, conf_t e5>
   struct mv<e1, e2, e3, e4, e5, -1, -1, -1, -1, -1, -1, -1>
   {
      typedef multivector<typename insert_element<e1,
                          typename insert_element<e2,
                          typename insert_element<e3,
                          typename insert_element<e4,
                          typename insert_element<e5,
                          cl_null>::clist>::clist>::clist>::clist>::clist, signature> type;
   };
   template<conf_t e1, conf_t e2, conf_t e3, conf_t e4>
   struct mv<e1, e2, e3, e4, -1, -1, -1, -1, -1, -1, -1, -1>
   {
      typedef multivector<typename insert_element<e1,
                          typename insert_element<e2,
                          typename insert_element<e3,
                          typename insert_element<e4,
                          cl_null>::clist>::clist>::clist>::clist, signature> type;
   };
   template<conf_t e1, conf_t e2, conf_t e3>
   struct mv<e1, e2, e3, -1, -1, -1, -1, -1, -1, -1, -1, -1>
   {
      typedef multivector<typename insert_element<e1,
                          typename insert_element<e2,
                          typename insert_element<e3,
                          cl_null>::clist>::clist>::clist, signature> type;
   };
   template<conf_t e1, conf_t e2>
   struct mv<e1, e2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1>
   {
      typedef multivector<typename insert_element<e1,
                          typename insert_element<e2,
                          cl_null>::clist>::clist, signature> type;
   };
   template<conf_t e1>
   struct mv<e1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1>
   {
      typedef multivector<typename insert_element<e1,
                          cl_null>::clist, signature> type;
   };
   /*template<>
   struct mv<-1, -1>
   {
      typedef multivector<cl_null, signature> type;
   };*/
};

} //end namespace gaalet

template<class A> inline
gaalet::multivector<typename A::clist, A::signature>
eval(const gaalet::expression<A>& a) {
   return gaalet::multivector<typename A::clist, A::signature>(a);
}


#endif
