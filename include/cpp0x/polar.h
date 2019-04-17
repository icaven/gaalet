#ifndef __GAALET_POLAR_H
#define __GAALET_POLAR_H

#include "utility.h"

namespace gaalet
{

template <class A>
struct polar : public expression<polar<A>> {
    static const conf_t I_conf = Power<2, A::metric::dimension>::value - 1;

    typedef configuration_list<I_conf, cl_null> clist;

    typedef typename A::metric metric;

    typedef typename A::element_t element_t;

    polar(const A& a_)
        : a(a_)
    {
    }

    template <conf_t... elements>
    auto element() const
    {
        static constexpr std::initializer_list<typename A::element_t> init = { 1 };
        return a * multivector<clist, metric, element_t>(init);
    }

protected:
    const A a;
};

} // end namespace gaalet

/// Polar of a multivector.
/** Implements the geometric product of an entity with the pseudoscalar.
 */
/// \ingroup ga_ops

template <class A>
inline gaalet::polar<A> polar(const gaalet::expression<A>& a)
{
    return gaalet::polar<A>(a);
}

#endif
