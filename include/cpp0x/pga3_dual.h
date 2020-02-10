#ifndef __GAALET_PGA3_DUAL_H
#define __GAALET_PGA3_DUAL_H
#pragma once

//
// Specialized version of the dual function for the PGA algebra
// Implements the Poincare duality
// See: section 2.3.1.3 of Gunn, Charles, "Geometry, Kinematics, and Rigid Body Mechanics in Cayley-Klein Geometries" [Gunn2011]
// The sign change is needed for some blades since the basis vectors in dual configuration
// are ordered to be in the canonical order instead of being permuted (as described by the referenced section)


#include "pga3.h"
#include "magnitude.h"
#include "utility.h"

namespace pga3 {

    namespace detail {
        template<conf_t I_conf, typename list, typename colist = cl_null>
        struct dual_list {
            typedef typename dual_list<I_conf, typename list::tail, typename insert_element<
                    I_conf ^ list::head, colist>::clist>::clist clist;
        };
        template<conf_t I_conf, typename colist>
        struct dual_list<I_conf, cl_null, colist> {
            typedef colist clist;
        };

        template<class A>
        struct dual : public expression<detail::dual<A >> {
            static const conf_t I_conf = Power<2, A::metric::dimension>::value - 1;

            typedef typename dual_list<I_conf, typename A::clist>::clist clist;

            typedef typename A::metric metric;

            typedef typename A::element_t element_t;

            dual(const A &a_)
                    : a(a_) {}

            template<conf_t conf>
            element_t element() const {
                // The odd permutation blades are stored negated, so change the sign when computing the dual
                // In the blade_conf() macro on the next lines, the order of the configuration list is always
                // shown in even parity order, but the actual ordering is either even or odd since the configuration is
                // just stored in a bitstring (ie. blade_conf() can't preserve the order)
                return (search_element<conf, clist>::index >= clist::size) ? 0.0 : a.template element<I_conf ^ conf>()
                                                                                   * ((conf == blade_conf(0, 2) || // even
                                                                                       conf == blade_conf(3, 1) || // odd
                                                                                       conf == vector_conf(1) || // even
                                                                                       conf == blade_conf(0, 3, 2) || // odd
                                                                                       conf == vector_conf(3) || // even
                                                                                       conf == blade_conf(0, 2, 1) // odd
                                                                                       ? -1 : 1)
                                                                                   );
            }

        protected:
            const A a;
        };

    } // end namespace detail

    /// @brief Dual of a PGA3 multivector.
    /// @ingroup ga_ops

    template <class A> inline
    detail::dual<A>
    dual(const gaalet::expression<A>& a) {
        return detail::dual<A>(a);
    }

} // end of namespace pga3

#endif // __GAALET_PGA3_DUAL_H