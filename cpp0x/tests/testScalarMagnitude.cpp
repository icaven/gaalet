#include "gaalet.h"

int main()
{
   gaalet::mv<1, 2, 4>::type a = {1, 2, 3};
   gaalet::mv<1, 2, 4>::type b = {3, 4, 5};

   gaalet::mv<0, 3, 5, 6>::type R = {cos(0.25*M_PI), sin(0.25*M_PI), 0.0, 0.0};


   std::cout << "scalar(a,a): " << scalar(a,a) << std::endl;
   
   std::cout << grade<0>(a*a) << std::endl;


   std::cout << "scalar(a*a, a): " << scalar(a*a, a) << std::endl;
}
