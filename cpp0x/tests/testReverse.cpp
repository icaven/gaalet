#include "gaalet.h"
#include <cmath>

int main()
{
   gaalet::mv<0, 3, 5, 6>::type R = {cos(-0.25*M_PI), sin(-0.25*M_PI), 0.0, 0.0};
   std::cout << "R: " << R << ", ~R: " << ~R << std::endl;

   gaalet::mv<1, 2, 4>::type a = { 1.0, 0.0, 0.0};

   std::cout << "a: " << a << ", R*a*(~R): " << grade<1>(R*a*(~R)) << std::endl;
}
