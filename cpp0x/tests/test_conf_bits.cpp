#include "gaalet.h"
#include "pga3.h"
#include <cmath>
#include <iostream>

using namespace gaalet;

// Projective Geometric Algebra 3d
// For a description of these operations, see:
// Gunn, Charles, “Projective geometric algebra: A new framework for doing euclidean geometry”,
// https://export.arxiv.org/abs/1901.05873
//
// Also, for a javascript implementation and description, see:
// https://github.com/enkimute/ganja.js

namespace test_conf_bits
{

template <int P, int Q, int R> 
struct space {
    typedef gaalet::algebra<signature<P, Q, R>, double> algebra;

    template <gaalet::conf_t... elements> struct mv {
        typedef typename algebra::template mv<elements...>::type type;
    };
};

// Convenience function and macro to generate the configuration list values
// Note that the null basis e0 has a bit representation that is expected by the gaalet implementation to square to 0
template <int P, int Q, int R>
inline constexpr conf_t vector_conf(const int v=-1) { 
    if (v < 0)
        return 0;
    else if (R > 0 && v == 0) {
        // Special case for the first degenerate base - allow it to be named with index 0
        return 1 << (P+Q+R-1);
    }
    else {
        return 1 << (v - 1); 
    }
};

template <int P, int Q, int R, int v_max=P+Q+R-(R>0?1:0)>
constexpr conf_t blade_conf(const int v=-1) {
    return vector_conf(v) | blade_conf<P, Q, R, v_max-1>(v);
};

template <int P, int Q, int R, int v_max>
constexpr conf_t blade_conf(const int v=-1) {
    return blade_conf<P, Q, R, v_max-1>(v);
};

template <int P, int Q, int R, 1>
constexpr conf_t blade_conf(const int v=-1) {
    return vector_conf(v);
}

#define blade_conf(v1, v2, ...) (pga3::vector_conf(v1) | pga3::vector_conf(v2) | pga3::vector_conf( __VA_ARGS__ ) )


// The basis elements of the dual algebra

// Scalar
const space::mv<0x0>::type one = { 1.0 };

// Vectors (representing planes)
const space::mv<vector_conf(0)>::type e0 = { 1.0 };
const space::mv<vector_conf(1)>::type e1 = { 1.0 };
const space::mv<vector_conf(2)>::type e2 = { 1.0 };
const space::mv<vector_conf(3)>::type e3 = { 1.0 };

// Bivectors (representing lines)
const space::mv<blade_conf(0, 1)>::type e01 = (e0 ^ e1);
const space::mv<blade_conf(0, 2)>::type e02 = (e0 ^ e2);
const space::mv<blade_conf(0, 3)>::type e03 = (e0 ^ e3);
const space::mv<blade_conf(1, 2)>::type e12 = (e1 ^ e2);
const space::mv<blade_conf(1, 3)>::type e13 = (e1 ^ e3);
const space::mv<blade_conf(2, 3)>::type e23 = (e2 ^ e3);
}
