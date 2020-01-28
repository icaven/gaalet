#ifndef __GAALET_STREAMING_H
#define __GAALET_STREAMING_H

#include <cmath>
#include "expression.h"
#include "multivector.h"

namespace gaalet {

        // These are specific to PGA3, and are ordered so that the configuration bit pattern can index into it
        // therefore they are out of the standard order
        static std::string basis_vector_names[16] = {"    ", "e1  ", "e2  ", "e12  ", "e3  ", "e13  ", "e23  ", "e123 ", "e0  ",
                                                     "e01  ", "e02  ", "e012 ", "e03  ", "e013 ", "e023 ", "e0123"};

//expression streaming
        template<typename G, typename clist>
        struct UnpackElementsToStream {
            template<class E, class T>
            //static void unpack(std::basic_ostream<E, T>& os, const G& e) {
            static void unpack(std::basic_ostream<E, T> &os, const gaalet::expression<G> &e_, bool previous_non_zero=false) {
                const G &e(e_);
                bool found_non_zero = previous_non_zero;
                if (e.template element<clist::head>() != 0) {
                    found_non_zero = true;
                    if (clist::head == 0) {
                        os << std::right << std::setw(8) << e.template element<0>();
                    } else {
//                        os << (std::signbit(e.template element<clist::head>()) ? "-" : " ");
//                        os << fabs(e.template element<clist::head>()) << basis_vector_names[clist::head];
                        os << std::right << std::setw(8) << e.template element<clist::head>() << basis_vector_names[clist::head];
                    }
                }

                UnpackElementsToStream<G, typename clist::tail>::unpack(os, e_, found_non_zero);
            }
        };


        template<typename G>
        struct UnpackElementsToStream<G, gaalet::cl_null> {
            template<class E, class T>
            static void unpack(std::basic_ostream<E, T> &os, const gaalet::expression<G> &, bool previous_non_zero=false)
            {
                if (!previous_non_zero) {
                    os << std::right << std::setw(8) << "0";
                }
            }
        };

        template<typename G, typename clist>
        struct UnpackConfigurationListToStream {
            template<class E, class T>
            static void unpack(std::basic_ostream<E, T> &os, const gaalet::expression<G> &e_, bool previous_non_zero=false) {
                const G &e(e_);
                bool non_zero_found = previous_non_zero;
                if (e.template element<clist::head>() != 0) {
                    non_zero_found = true;
                    os << std::hex << clist::head << " ";
                }
                UnpackConfigurationListToStream<G, typename clist::tail>::unpack(os, e_, non_zero_found);
            }
        };

        template<typename G>
        struct UnpackConfigurationListToStream<G, gaalet::cl_null> {
            template<class E, class T>
            static void unpack(std::basic_ostream<E, T> &os, const gaalet::expression<G> &, bool previous_non_zero=false)
            {
                if (!previous_non_zero) {
                    os << "0 ";
                }
            }
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

//   os << "[ " << std::dec;
      gaalet::UnpackElementsToStream<G, typename G::clist>::unpack(os, e);
//   os << "] ";
//   os << " { " << std::hex;
//      gaalet::UnpackConfigurationListToStream<G, typename G::clist>::unpack(os, e);
//   os << "}";

   return os;
}

template<class E, class T, class G>
std::basic_ostream<E, T>& operator<<(std::basic_ostream<E, T>&& os, const gaalet::expression<G>& e)
{
   return (os << e);
}

}  // end namespace std

#endif
