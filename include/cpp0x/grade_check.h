#ifndef __GAALET_GRADE_CHECK_H
#define __GAALET_GRADE_CHECK_H

#pragma once

#include "configuration_list.h"
#include "utility.h"

namespace gaalet {
    // check for configuration list that describes a scalar only
    template<typename CL>
    struct check_scalar {
        static const bool value = CL::size == 1 && (BitCount<CL::head>::value == 0);
    };

    // check for configuration list with vectors only
    template<typename CL>
    struct check_vector {
        static const bool value = (BitCount<CL::head>::value == 1) ? check_vector<typename CL::tail>::value : false;
    };
    template<>
    struct check_vector<cl_null> {
        static const bool value = true;
    };

    // check for configuration list with bivectors only
    template<typename CL>
    struct check_bivector {
        static const bool value = (BitCount<CL::head>::value == 2) ? check_bivector<typename CL::tail>::value : false;
    };
    template<>
    struct check_bivector<cl_null> {
        static const bool value = true;
    };

    // check for configuration list with trivectors only
    template<typename CL>
    struct check_trivector {
        static const bool value = (BitCount<CL::head>::value == 3) ? check_trivector<typename CL::tail>::value : false;
    };
    template<>
    struct check_trivector<cl_null> {
        static const bool value = true;
    };

    // check for configuration list that describes a 4-vectors only
    template<typename CL>
    struct check_quadvector {
        static const bool value = (BitCount<CL::head>::value == 4) ? check_quadvector<typename CL::tail>::value : false;
    };
    template<>
    struct check_quadvector<cl_null> {
        static const bool value = true;
    };

    //check that all grades in the multivector are even
    template<typename CL>
    struct check_even_grade
    {
        static const bool value = (BitCount<CL::head>::value % 2 == 0) ? check_even_grade<typename CL::tail>::value : false;
    };
    template<>
    struct check_even_grade<cl_null>
    {
        static const bool value = true;
    };

    // A blade is a outer product of vectors, so if the configuration list is of length 1
    // and if the bitcount of entry isn't 0 (which would indicate a scalar), then the entity is a blade.
    template<typename CL>
    struct check_blade {
        static const bool value = CL::size == 1 && BitCount<CL::head>::value > 0;
    };

} // end of namespace gaalet

#endif // __GAALET_GRADE_CHECK_H
