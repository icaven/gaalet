#include "gaalet.h"
#include <sys/time.h>

int main()
{
   timeval start, end;

   gaalet::fast_mv<0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16>::type a = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};
   gaalet::fast_mv<1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16>::type b = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};

   gaalet::fast_mv<2, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16>::type c;

   gettimeofday(&start, 0);
   for(int i = 0; i<1e7; ++i) {
      c = c + a + b - a - b + a + b - a - b - c;
   }
   gettimeofday(&end, 0);
   double solveTime = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)*1e-6;

   //std::cout << "a: " << a << std::endl;
   //std::cout << "b: " << b << std::endl;
   //std::cout << "c: " << c << std::endl;

   std::cout << "add solve time: " << solveTime << std::endl;

   gaalet::mv<3, 1, 2, 4>::type d = {1.0, 2.0, 3.0};
   gaalet::mv<4, 1, 2, 4>::type e = {5.0, 6.0, 7.0};
   gaalet::mv<5, 0, 3, 5, 6>::type f = {8.0, 9.0, 10.0, 11.0};

   gettimeofday(&start, 0);
   for(int i = 0; i<1e7; ++i) {
      f = d * e * f;
   }
   gettimeofday(&end, 0);
   solveTime = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)*1e-6;

   //std::cout << "d: " << d << std::endl;
   //std::cout << "e: " << e << std::endl;
   //std::cout << "f: " << f << std::endl;
   std::cout << "multiply solve time: " << solveTime << std::endl;
}
