#include "gaalet.h"

typedef gaalet::algebra<gaalet::signature<0,0> > em;

int main()
{
   em::mv<3, 5>::type a(1.0, 2.0);
   em::mv<3, 4>::type b(4.0, 7.0);

   typeof(a-b) c = a - b;

   std::cout << "a: " << a << std::endl;
   std::cout << "b: " << b << std::endl;
   std::cout << "c: " << c << std::endl;

   em::mv<0>::type d;
}
