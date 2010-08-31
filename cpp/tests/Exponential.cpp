#define _USE_MATH_DEFINES
#include "gaalet.h"


typedef gaalet::algebra<gaalet::signature<3,0> > em;
typedef gaalet::algebra<gaalet::signature<4,1> > cm;

int main()
{
   em::mv<3, 5, 6>::type A(3, 4, 5);
   std::cout << "A: " << A << ", exp(A): " << exp(A) << std::endl; 
  
   em::mv<1, 2, 4>::type a(1.0, 2.0, 3.0);
   em::mv<3, 5, 6>::type m(-0.25*M_PI, 0.0, 0.0);
   em::mv<0, 3, 5, 6>::type R = exp(m);
   std::cout << "R: " << m << ", R=exp(m): " << R << std::endl;
   std::cout << "a: " << a << ", R*a*(~R): " << R*a*(~R) << ", R*a*(!R):" << R*a*(!R) << std::endl;

   cm::mv<0x01>::type e1(1.0);
   cm::mv<0x02>::type e2(1.0);
   cm::mv<0x04>::type e3(1.0);
   cm::mv<0x08>::type ep(1.0);
   cm::mv<0x10>::type em(1.0);

   cm::mv<0x08, 0x10>::type e0 = 0.5*(em-ep);
   cm::mv<0x08, 0x10>::type einf = em+ep;

   cm::mv<0x08^1,0x08^2,0x08^4,0x10^1,0x10^2,0x10^4>::type S = einf*(2.0*e1 + 1.0*e2 + 0.5*e3);
   std::cout << "S: " << S << ", S*S: " << S*S << std::endl;
   cm::mv<0x00,0x08^1,0x08^2,0x08^4,0x10^1,0x10^2,0x10^4>::type T = exp(0.5*S);
   std::cout << "T: " << T << ", T*e0*(~T): " << T*e0*(~T) << std::endl;
   std::cout << "T: " << T << ", <T*e0*(~T)>_1: " << grade<1>(T*e0*(~T)) << std::endl;

   //vvv fails, because a is no bivector
   //std::cout << exp(a) << std::endl;
}
