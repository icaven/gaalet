#ifndef __GAALET_PGA3_OPS_H
#define __GAALET_PGA3_OPS_H

#pragma once

#include "gaalet.h"
#include "pga3.h"
#include "pga3_dual.h"
#include "pga3_normalize.h"
#include "pga3_point.h"
#include "pga3_line.h"
#include "pga3_plane.h"

//
// PGA3 geometric operations
//

namespace pga3 {

    template<typename X>
    inline
    auto polar(const gaalet::expression<X> &x) {
    //    return x * I;
        // Need to add the zero scalar so that the configuration list won't be empty after some geometric products
        return 0.0 * one + x * I;
    }

    // Sandwich product
    template<typename L, typename R>
    inline
    auto sandwich(const gaalet::expression<L> &l, const gaalet::expression<R> &r) {
        return (r * l * (~r));
    }

    // Regressive product
    template<typename L, typename R>
    inline
    auto vee(const gaalet::expression<L> &l, const gaalet::expression<R> &r) {
        return pga3::dual(pga3::dual(l) ^ pga3::dual(r));
    }

    // Commutator product
    template<typename L, typename R>
    inline
    auto commutator(const gaalet::expression<L> &l, const gaalet::expression<R> &r) {
        return 0.5 * (l * r - r * l);
    }

    // Anti-commutator product
    template<typename L, typename R>
    inline
    auto anti_commutator(const gaalet::expression<L> &l, const gaalet::expression<R> &r) {
        return 0.5 * (l * r + r * l);
    }

    // Joins
    template<typename L, typename R>
    inline
    auto line_from_points(const gaalet::expression<L> &start, const gaalet::expression<R> &end) {
        return vee(start, end);
    }

    template<class T1, class T2, class T3>
    inline
    auto plane_from_points(const T1 &P1, const T2 &P2, const T3 &P3) {
        return pga3::normalize(pga3::vee(pga3::vee(P1, P2), P3));
    }

    template<class TP, class TL>
    inline
    auto plane_from_point_and_line(const TP &P, const TL &L) {
        return vee(P, L);
    }

    // Meets
    template<class T1, class T2>
    inline
    auto line_from_planes(const T1 &p1, const T2 &p2) {
        return p1 ^ p2;
    }

    template<class T1, class T2, class T3>
    inline
    auto point_from_planes(const T1 &p1, const T2 &p2, const T3 &p3) {
        return pga3::normalize(p1 ^ p2 ^ p3);
    }

    template<class TP, class TL>
    inline
    auto point_from_line_and_plane(const TL &L, const TP &P) {
        return pga3::normalize(L ^ P);
    }

    // Points are in a CCW direction
    template<class T1, class T2, class T3>
    inline
    auto normal_to_plane(const T1 &P1, const T2 &P2, const T3 &P3) {
        return pga3::normalize(P1 & plane_from_points(P1, P2, P3));
    }

    //
    // Rigid motions
    //

    // Rotation around an axis
    template<class A, class T>
    inline
    auto rotor(const A &axis, const T &angle) {
        return one * cos(-angle * 0.5) + sin(-angle * 0.5) * pga3::normalize(axis);
    };

    // Translation in the direction of a line
    template<class L, class T>
    inline
    auto translator(const L &line, const T &distance) {
        return one + 0.5 * distance * pga3::normalize(line) * I;
    }

    // Motor operation implements a rotation followed by a translation, which is similar to
    // a screw operation, which is a rotation around a line while travelling a distance along the line
    template<class L, class T, typename P = space::algebra::element_t>
    inline
    auto motor(const L &line, const T &distance, const P angle) {
        return translator(line, distance) * rotor(line, angle);
    }

    //template <class L, class T, typename P = space::algebra::element_t>
    //auto screw(const L& line, const T& distance, const P pitch)
    //{
    //    return exp(0.5 * distance * (pga3::normalize(line) + pitch * pga3::normalize(line) * I));
    //}

    // Alternative form of the screw operation, to be used for comparison
    template<class L, class T, typename P = space::algebra::element_t>
    auto screw(const L &line, const T &distance, const P pitch) {
        return 0.5 * distance * (one + pitch * I) * pga3::normalize(line);
    }

} // end of namespace pga3

#endif // __GAALET_PGA3_OPS_H
