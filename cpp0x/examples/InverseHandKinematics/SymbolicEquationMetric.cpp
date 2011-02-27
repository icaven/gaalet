#include "gaalet.h"

int main() {
   typedef gaalet::algebra< gaalet::signature<3,0>, gaalet::symbex> sem;
   sem::mv<0>::type one = {1.0};

   sem::mv<1,2,4>::type x = {"x1", "x2", "x3"};
   sem::mv<1,2,4>::type t = {"t1", "t2", "t3"};

   sem::mv<3,5,6>::type m = {"m12", "m13", "m23"};
   auto mag_m = magnitude(m);
   auto R = eval(cos(eval((-0.5)*mag_m))*one + m*!mag_m*sin(eval((-0.5)*mag_m)));
   //std::cout << "R*~R: " << grade<0>(R*~R) << std::endl;
   std::cout << "t&(R*x*~R): " << grade<0>(t&(R*x*~R)) << std::endl;

   auto U = eval(0.5*(t - grade<1>(R*x*~R))&(t - grade<1>(R*x*~R)));
   //std::cout << "U: " << U << std::endl;
}
