#include "gaalet.h"

typedef gaalet::algebra<gaalet::signature<0,3,1>> ma;

int main()
{
   ma::mv<1>::type e1={1.0};
   ma::mv<2>::type e2={1.0};
   ma::mv<4>::type e3={1.0};
   ma::mv<8>::type e0={1.0};


   auto x = ma::mv<0>::type({1.0}) + 1.0*e1*e0 + 2.0*e2*e0 + 3.0*e3*e0;
   std::cout << "x: " << x << std::endl;

   double phi = M_PI*0.5;
   auto r = ma::mv<0>::type({cos(phi*0.5)}) + sin(phi*0.5)*(0.0*e2*e3 + 0.0*e3*e1 + 1.0*e1*e2);
   std::cout << "r: " << r << std::endl;

   auto t = ma::mv<0>::type({1.0}) + 0.5*(10.0*e1*e0 + 0.0*e2*e0 + 0.0*e3*e0);
   std::cout << "t: " << t << std::endl;

   std::cout << "r*x*~r: " << r*x*~r << std::endl;
   std::cout << "t*x*~t: " << t*x*~t << std::endl;

   auto g = t*r;
   std::cout << "g: " << g << std::endl;

   std::cout << "g*x*~g: " << g*x*~g << std::endl;
}
