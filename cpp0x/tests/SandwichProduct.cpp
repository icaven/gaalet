#include "gaalet.h"

template<typename L, typename R> inline
auto sandwich(const gaalet::expression<L>& l_, const gaalet::expression<R>& r_) -> decltype(R()*L()*(~R()))
{
   const L& l(l_);
   const R& r(r_);

   return (r*l*(~r));
}


typedef gaalet::algebra< gaalet::signature<3,0> > em;

int main()
{
   em::mv<1,2,4>::type v = {1.0, 2.0, 3.0};
   em::mv<0,3,5,6>::type R = {cos(-0.5*0.5*M_PI), sin(-0.5*0.5*M_PI), 0.0, 0.0};

   std::cout << "R*v*~R: " << R*v*~R << std::endl;
   std::cout << "<R*v*~R>_1: " << grade<1>(R*v*~R) << std::endl;

   std::cout << "sandwich(v,R): " << sandwich(v,R) << std::endl;
}
