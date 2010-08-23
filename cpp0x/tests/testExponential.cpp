#include "gaalet.h"

typedef gaalet::algebra<gaalet::signature<3,0>> em;

int main()
{
   em::mv<3, 5, 6>::type A = {3, 4, 5};
   std::cout << "A: " << A << ", exp(A): " << exp(A) << std::endl; 
  
   em::mv<1, 2, 4>::type a = {1.0, 2.0, 3.0};
   em::mv<3, 5, 6>::type m = {-0.25*M_PI, 0.0, 0.0};
   auto R = exp(m);
   std::cout << "R: " << m << ", R=exp(m): " << R << std::endl;
   std::cout << "a: " << a << ", R*a*(~R): " << R*a*(~R) << ", R*a*(!R):" << R*a*(!R) << std::endl;
}
