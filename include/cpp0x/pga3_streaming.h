#ifndef __GAALET_PGA3_STREAMING_H
#define __GAALET_PGA3_STREAMING_H

#pragma once

#include <ostream>
#include <iostream>

#include "gaalet.h"
#include "pga3.h"

namespace pga3 {

// These functions may be used to print out PGA3 multivectors in a more descriptive way than just
// using the functions in streaming.h.
// @todo: Find a way to call these functions for gaalet::expressions that contain pga3 multivectors


// These are specific to PGA3, and are ordered so that the configuration bit pattern can index into it
static std::string basis_vector_names[16] = {"    ", "e0  ", "e1  ", "e01 ", "e2  ", "e02 ", "e12 ", "e021 ", "e3  ",
                                             "e03 ", "e31 ", "e013 ", "e23 ", "e032 ", "e123 ", "e0123"};

    // PGA3 multivector streaming

    template <typename clist>
    struct UnpackElementsToStream {
        template<class E, class T>
        static void unpack(std::basic_ostream<E, T> &os, const gaalet::multivector<clist,
                pga3::space::algebra::metric, pga3::space::algebra::element_t>& e,
                bool previous_non_zero=false) {
            bool found_non_zero = previous_non_zero;
            if (e.template element<clist::head>() != 0) {
                found_non_zero = true;
                if (clist::head == 0) {
                    os << std::right << std::setw(8) << e.template element<0>();
                }
                else if (clist::head == J_CONF || clist::head == E1_CONF || clist::head == E3_CONF) {
                    // The odd permutation blades are stored with negated values, so change the sign
                    os << std::right << std::setw(8) << -1 * e.template element<clist::head>() << basis_vector_names[clist::head];
                }
                else {
                    os << std::right << std::setw(8) << e.template element<clist::head>() << basis_vector_names[clist::head];
                }
            }
            UnpackElementsToStream<typename clist::tail>::unpack(os, e, found_non_zero);
        }

    };

    template<>
    struct UnpackElementsToStream<gaalet::cl_null> {
        template<class E, class T>
        static void unpack(std::basic_ostream<E, T> &os, const gaalet::multivector<gaalet::cl_null,
                pga3::space::algebra::metric, pga3::space::algebra::element_t>& , bool previous_non_zero=false)
        {
            if (!previous_non_zero) {
                os << std::right << std::setw(8) << "0";
            }
        }
    };

} // end of namespace pga3

namespace std {

    template<typename CL>
    std::ostream &operator<<(std::ostream &os, const gaalet::multivector<CL,
            pga3::space::algebra::metric, pga3::space::algebra::element_t> &e) {
        pga3::UnpackElementsToStream<CL>::unpack(os, e);
        return os;
    }

    template<class E, class T, class CL>
    std::basic_ostream<E, T> &operator<<(std::basic_ostream<E, T> &&os,
                                         const gaalet::expression<gaalet::multivector<CL, pga3::space::algebra::metric,
                                                 pga3::space::algebra::element_t>> &e) {
        pga3::UnpackElementsToStream<CL>::unpack(os, e);
        return os;
    }

} // end of namespace std

#endif // __GAALET_PGA3_STREAMING_H
