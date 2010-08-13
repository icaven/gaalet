#include "gaalet.h"

int main()
{
   gaalet::mv<1, 2, 4>::type a = {1, 2, 3};

   //std::cout << "a: " << a << ", !a: " << !a << std::endl;
   std::cout << "a: " << a << std::endl;
   std::cout << "!a: " << (~a)*(1.0/(a*(~a)).element<0x00>()) << std::endl;
   std::cout << "!a: " << (!a) << std::endl;
}
