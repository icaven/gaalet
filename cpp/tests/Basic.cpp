#include "gaalet.h"


int main()
{
   gaalet::mv<3, 5>::type a = {1.0, 2.0};
   gaalet::mv<3, 4>::type b = {4.0, 7.0};

   auto c = a - b;

   std::cout << "a: " << a << std::endl;
   std::cout << "b: " << b << std::endl;
   std::cout << "c: " << c << std::endl;

   gaalet::mv<0>::type d;
}
