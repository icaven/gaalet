#ifndef __GAALET_PGA3_H
#define __GAALET_PGA3_H

#pragma once

// Projective Geometric Algebra 3d
// For a description of these operations, see [Gunn2019]:
// Gunn, Charles, “Projective geometric algebra: A new framework for doing euclidean geometry”,
// https://export.arxiv.org/abs/1901.05873
//
// and:
// Gunn, Charles, "Geometry, Kinematics, and Rigid Body Mechanics in Cayley-Klein Geometries" [Gunn2011]


//
// Also, for a javascript implementation and description, see:
// https://github.com/enkimute/ganja.js

#include "gaalet.h"

using namespace gaalet;

namespace pga3 {
    ///@brief Definition of the 4 dimensional geometric algebra space that is projective
    ///       and uses dual basis elements and has a multivector type of variable configuration
    struct space {
        typedef gaalet::algebra<signature<3, 0, 1>, double> algebra;

        template<gaalet::conf_t... elements>
        struct mv {
            typedef typename algebra::mv<elements...>::type type;
        };
    };

    // Convenience function and macro to generate the configuration list values
    // Note that the order of the vectors in the blade_conf is not retained (because the result is just a bit string),
    // but for consistency they will be shown in the order of the basis vectors in the blades
    inline constexpr conf_t vector_conf(const int v = -1) { return v < 0 ? 0 : 1 << v; };
    #define blade_conf(v1, v2, ...) (pga3::vector_conf(v1) | pga3::vector_conf(v2) | pga3::vector_conf( __VA_ARGS__ ) )


    // The basis elements of the dual algebra

    // Scalar
    const space::mv<0x0>::type one = {1.0};

    // Vectors (representing planes)
    const space::mv<vector_conf(0)>::type e0 = {1.0};
    const space::mv<vector_conf(1)>::type e1 = {1.0};
    const space::mv<vector_conf(2)>::type e2 = {1.0};
    const space::mv<vector_conf(3)>::type e3 = {1.0};

    // Bivectors (representing lines)
    const space::mv<blade_conf(0, 1)>::type e01 = (e0 ^ e1);
    const space::mv<blade_conf(0, 2)>::type e02 = (e0 ^ e2);
    const space::mv<blade_conf(0, 3)>::type e03 = (e0 ^ e3);
    const space::mv<blade_conf(1, 2)>::type e12 = (e1 ^ e2);
    const space::mv<blade_conf(1, 3)>::type e31 = (e3 ^ e1);
    const space::mv<blade_conf(2, 3)>::type e23 = (e2 ^ e3);

    // Alternative names for biquaternion components
    constexpr conf_t DUAL_I_CONF = blade_conf(0, 3);
    constexpr conf_t DUAL_J_CONF = blade_conf(0, 2);
    constexpr conf_t DUAL_K_CONF = blade_conf(0, 1);
    constexpr conf_t I_CONF = blade_conf(1, 2);
    constexpr conf_t J_CONF = blade_conf(3, 1);
    constexpr conf_t K_CONF = blade_conf(2, 3);
    const space::mv<I_CONF>::type i = e12;      // z-axis      - see Section 7.4.3 of [Gunn2011]
    const space::mv<J_CONF>::type j = e31;      // y-axis
    const space::mv<K_CONF>::type k = e23;      // x-axis
    const space::mv<DUAL_K_CONF>::type dual_k = e01;
    const space::mv<DUAL_J_CONF>::type dual_j = e02;
    const space::mv<DUAL_I_CONF>::type dual_i = e03;


    // The tri-vectors (representing points)
    constexpr conf_t E0_CONF = blade_conf(1, 2, 3);
    constexpr conf_t E1_CONF = blade_conf(0, 3, 2);
    constexpr conf_t E2_CONF = blade_conf(0, 1, 3);
    constexpr conf_t E3_CONF = blade_conf(0, 2, 1);


    const space::mv<E0_CONF>::type e123 = (e1 ^ e2 ^ e3);
    const space::mv<E1_CONF>::type e032 = (e0 ^ e3 ^ e2);
    const space::mv<E2_CONF>::type e013 = (e0 ^ e1 ^ e3);
    const space::mv<E3_CONF>::type e021 = (e0 ^ e2 ^ e1);

    // Alternative names
    const space::mv<E0_CONF>::type E0 = e123; // This is the origin point in Euclidean space
    const space::mv<E1_CONF>::type E1 = e032;
    const space::mv<E2_CONF>::type E2 = e013;
    const space::mv<E3_CONF>::type E3 = e021;

    // The pseudoscalar (representing all space)
    const gaalet::conf_t pseudoscalar_conf = (1 << space::algebra::metric::dimension) - 1;
    const space::mv<pseudoscalar_conf>::type I = (e0 ^ e1 ^ e2 ^ e3);

    // Types for common geometric entities
    typedef space::mv<E0_CONF, E1_CONF, E2_CONF, E3_CONF>::type Point_t;
    typedef space::mv<DUAL_K_CONF, DUAL_J_CONF, DUAL_I_CONF, I_CONF, J_CONF, K_CONF>::type Line_t;
    typedef space::mv<vector_conf(1), vector_conf(2), vector_conf(3), vector_conf(0)>::type Plane_t;
    typedef space::mv<0, DUAL_K_CONF, DUAL_J_CONF, DUAL_I_CONF, I_CONF, J_CONF, K_CONF, pseudoscalar_conf>::type Motor_t;

} // end namespace pga3

#endif // __GAALET_PGA3_H
