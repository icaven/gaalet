#include "gaalet.h"
#include <iostream>
#include <cstdlib>

typedef gaalet::algebra<gaalet::signature<4,1> > cm;

__global__ void test()
{
   double r = (double)blockDim.x;
   
   unsigned int n_q = blockDim.x*blockDim.y*blockDim.z;
   __shared__ unsigned int n_s;
   if(threadIdx.x == 0 && threadIdx.y==0 && threadIdx.z==0) n_s = 0;

   cm::mv<0x01>::type e1(1.0);
   cm::mv<0x02>::type e2(1.0);
   cm::mv<0x04>::type e3(1.0);
   cm::mv<0x08>::type ep(1.0);
   cm::mv<0x10>::type em(1.0);

   cm::mv<0x00>::type one(1.0);

   cm::mv<0x08, 0x10>::type e0 = 0.5*(em-ep);
   cm::mv<0x08, 0x10>::type einf = em+ep;

   cm::mv<0x08, 0x10>::type S = eval(e0 - 0.5*r*r*einf);


   cm::mv<0x01, 0x02, 0x04>::type x = eval(((double)threadIdx.x*e1 + (double)threadIdx.y*e2 + (double)threadIdx.z*e3)*r);
   cm::mv<0x01, 0x02, 0x04, 0x08, 0x10>::type P = x + 0.5*(x&x)*einf + e0;
   double d = eval(S&P);
   if(d>=0.0) {
      atomicAdd(&n_s, 1);
   }
   __syncthreads();
}


int main()
{
   std::cout << "Hello Gaalet Monte Carlo on Cuda!" << std::endl;

   cudaSetDevice( 0 );

   dim3 threads( 10, 10, 10 );

   test <<< 1, threads >>>();

   cudaThreadExit();

   //std::cout << "Pi: " << 6.0*(double)n_s/(double)n_q << std::endl;
}
