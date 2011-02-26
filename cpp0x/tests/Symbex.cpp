#include "gaalet.h"

int main()
{
   gaalet::symbex a("a");
   gaalet::symbex b("b");

   gaalet::symbex c = a + b;

   std::cout << "a: " << a << ", b: " << b << ", c: " << c << ", a+c: " << a+c << ", (a+c)*b: " << (a+c)*b << std::endl;


   typedef gaalet::algebra< gaalet::signature<0,0>, gaalet::symbex> sr;

   sr::mv<0>::type sa = {"a"};
   sr::mv<0>::type sb = {"b"};

   auto sc = a + b;
   
   std::cout << "a: " << sa << ", b: " << sb << ", c: " << sc << ", a+c: " << sa+sc << ", (a+c)*b: " << (sa+sc)*sb << ", sa.element<1>(): " << sa.element<1>() << std::endl;


   typedef gaalet::algebra< gaalet::signature<3,0>, gaalet::symbex> sem;
   sem::mv<1,2,4>::type x = {"x1", "x2", "x3"};
   sem::mv<0,3,5,6>::type R = {"R0", "R12", "R13", "R23"};

   std::cout << "x: " << x << ", R: " << R << ", R*x*~R: " << R*x*~R << std::endl;
}
