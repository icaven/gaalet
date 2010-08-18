#include "gaalet.h"
#include <sys/time.h>
#include <cmath>

int main()
{
   timeval start, end;
   double solveTime;

   typedef gaalet::algebra<gaalet::signature<12,0> > em;

   em::mv<1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12>::type a;
   em::mv<1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12>::type b;
   
   em::mv<1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12>::type c;

   gettimeofday(&start, 0);
   for(int i = 0; i<1e7; ++i) {
      c = c + a + b - a - b + a + b - a - b - c;
   }
   gettimeofday(&end, 0);
   solveTime = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)*1e-6;

   std::cout << "a: " << a << std::endl;
   std::cout << "b: " << b << std::endl;
   std::cout << "c: " << c << std::endl;

   std::cout << "operator=(): add solve time: " << solveTime << std::endl;


   gettimeofday(&start, 0);
   for(int i = 0; i<1e7; ++i) {
      c = eval(c + a + b - a - b + a + b - a - b - c);
   }
   gettimeofday(&end, 0);
   solveTime = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)*1e-6;

   std::cout << "a: " << a << std::endl;
   std::cout << "b: " << b << std::endl;
   std::cout << "c: " << c << std::endl;

   std::cout << "eval(): add solve time: " << solveTime << std::endl;


   gettimeofday(&start, 0);
   for(int i = 0; i<1e7; ++i) {
      c.assign(c + a + b - a - b + a + b - a - b - c);
   }
   gettimeofday(&end, 0);
   solveTime = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)*1e-6;

   std::cout << "a: " << a << std::endl;
   std::cout << "b: " << b << std::endl;
   std::cout << "c: " << c << std::endl;

   std::cout << "assign(): add solve time: " << solveTime << std::endl;
}
