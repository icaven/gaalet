#ifndef __GAALET_PGA3_UTILITY_H
#define __GAALET_PGA3_UTILITY_H

#pragma once

#include <ostream>
#include <cmath>
#include <cfloat>

#include "gaalet.h"
#include "pga3.h"

namespace pga3 {

    inline bool isclose(double a, double b, double rtol = 1e-05, double atol = DBL_EPSILON) {
        return fabs(a - b) <= (atol + rtol * fabs(b));
    }

    /// Return the Study number and the axis (a normalized bivector) associated with the given bivector
    // See: Section 7.7 of Gunn, Charles, "Geometry, Kinematics, and Rigid Body Mechanics in Cayley-Klein Geometries" [Gunn2011]
    // Normalize the bivector by multiplying it by the inverse of the square root of the associated Study number z, and the result
    // is B, which is the axis of the given bivector
    std::pair<std::pair<double, double>, Line_t> inline
    bivector_axis(pga3::Line_t bivector) {
        auto z = ::eval(bivector * (~bivector));
        double a = ::grade<0>(z).template element<0>();
        auto axis = ::eval(
                bivector);   // In the case that the inverse sqrt of z doesn't exist (when a == 0), the B is equal to the bivector
        double b = 0.0;
        std::pair<double, double> study_number = {0.0, 0.0};
        if (!isclose(a, 0.0)) {
            auto sqrt_a = sqrt(a);
            b = ::grade<4>(z).template element<pseudoscalar_conf>();
            auto inv_sqrt_z = sqrt_a / a * one - ((b) / (2 * a * sqrt_a)) * I;
            axis = inv_sqrt_z * axis;
            study_number = std::make_pair(sqrt_a, b / (2.0 * sqrt_a));
            return std::make_pair(study_number, axis);
        }
        return std::make_pair(study_number, axis);
    }

} // end of namespace pga3

#endif // __GAALET_PGA3_UTILITY_H
