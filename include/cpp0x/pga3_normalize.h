#ifndef __GAALET_PGA3_NORMALIZE_H
#define __GAALET_PGA3_NORMALIZE_H

#pragma once

// Normalize a PGA3 multivector

#include "pga3_norm.h"

namespace pga3 {

    template<class A>
    inline
    gaalet::scalar_multivector_product<A>
    normalize(const gaalet::expression<A> &a) {
        return 1.0/pga3::norm(a) * a;
    }

} // end of namespace pga3

#endif // __GAALET_PGA3_NORMALIZE_H