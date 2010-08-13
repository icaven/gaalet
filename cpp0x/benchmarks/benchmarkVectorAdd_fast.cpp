#include "gaalet.h"
#include <sys/time.h>
#include <cmath>

int main()
{
   timeval start, end;
   double solveTime;

   gaalet::fast_mv<0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16>::type a = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};
   gaalet::fast_mv<1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16>::type b = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};

   gaalet::fast_mv<2, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16>::type c;

   gettimeofday(&start, 0);
   for(int i = 0; i<1e7; ++i) {
      c = c + a + b - a - b + a + b - a - b - c;
   }
   gettimeofday(&end, 0);
   solveTime = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)*1e-6;

   //std::cout << "a: " << a << std::endl;
   //std::cout << "b: " << b << std::endl;
   //std::cout << "c: " << c << std::endl;

   std::cout << "fast_mv: operator=(): add solve time: " << solveTime << std::endl;


   gaalet::fast<3>::mv<1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16>::type d = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};
   gaalet::fast<4>::mv<1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16>::type e = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};

   gaalet::fast<5>::mv<1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16>::type f;

   gettimeofday(&start, 0);
   for(int i = 0; i<1e7; ++i) {
      f = f + d + e - d - e + d + e - d - e - f;
   }
   gettimeofday(&end, 0);
   solveTime = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)*1e-6;

   std::cout << "mv::fast: operator=(): add solve time: " << solveTime << std::endl;
}
