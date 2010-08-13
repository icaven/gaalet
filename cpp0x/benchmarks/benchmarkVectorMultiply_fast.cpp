#include "gaalet.h"
#include <sys/time.h>
#include <cmath>

int main()
{
   timeval start, end;

   gaalet::fast_mv<0, 1, 2, 4, 7>::type a = {1.0, 0.0, 0.0, 0.0};
   gaalet::fast_mv<1, 0, 3, 5, 6>::type b = {cos(-M_PI*0.25), sin(-M_PI*0.25), 0.0, 0.0};

   gettimeofday(&start, 0);
   for(int i = 0; i<1e8; ++i) {
      //a = (~b)*b*a*(~b)*b;
      a = b*a*(~b);
   }
   gettimeofday(&end, 0);
   double solveTime = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)*1e-6;

   //std::cout << "d: " << d << std::endl;
   //std::cout << "e: " << e << std::endl;
   //std::cout << "f: " << f << std::endl;
   std::cout << "a: " << "( " << a[0] << " " << a[1] << " " << a[2] << " " << a[3] << ") " << std::endl;
   std::cout << "b: " << "( " << b[0] << " " << b[1] << " " << b[2] << " " << b[3] << ") " << std::endl;
   std::cout << "fast_mv: operator=(): multiply solve time: " << solveTime << std::endl;


   gaalet::fast<2>::mv<1, 2, 4, 7>::type d = {1.0, 0.0, 0.0, 0.0};
   gaalet::fast<3>::mv<0, 3, 5, 6>::type e = {cos(-M_PI*0.25), sin(-M_PI*0.25), 0.0, 0.0};

   gettimeofday(&start, 0);
   for(int i = 0; i<1e8; ++i) {
      //d = (~e)*e*d*(~e)*e;
      d = e*d*(~e);
   }
   gettimeofday(&end, 0);
   solveTime = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)*1e-6;

   //std::cout << "d: " << d << std::endl;
   //std::cout << "e: " << e << std::endl;
   //std::cout << "f: " << f << std::endl;
   std::cout << "d: " << "( " << d[0] << " " << d[1] << " " << d[2] << " " << d[3] << ") " << std::endl;
   std::cout << "e: " << "( " << e[0] << " " << e[1] << " " << e[2] << " " << e[3] << ") " << std::endl;
   std::cout << "mv::fast: operator=(): multiply solve time: " << solveTime << std::endl;
}
