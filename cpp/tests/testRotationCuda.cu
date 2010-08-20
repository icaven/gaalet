#include "gaalet.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <cmath>

typedef gaalet::algebra<gaalet::signature<3,0> > em;
typedef em::mv<0x01,0x02,0x04>::type Vector;
typedef em::mv<0x00,0x03,0x05,0x06>::type Rotor;

struct rotation_functor
{
   rotation_functor(const Rotor& setR)
      : R(setR),
        invR(!setR)
   { }

   __host__ __device__
      Vector operator()(const Vector& x) const
      { 
         return grade<1>(R*x*invR);
      }

   Rotor R;
   Rotor invR;
};


int main()
{
   thrust::host_vector<Vector> h_x(100);
   h_x[0][0] = 1.0; h_x[0][1] = 0.0; h_x[0][2] = 0.0;
   h_x[1][0] = 0.0; h_x[1][1] = 1.0; h_x[1][2] = 0.0;
   h_x[2][0] = 0.0; h_x[2][1] = 0.0; h_x[2][2] = 1.0;

   thrust::device_vector<Vector> d_x = h_x;

   thrust::device_vector<Vector> d_y(100);

   Rotor R;
   R[0] = cos(-0.5*0.5*M_PI); R[1] = sin(-0.5*0.5*M_PI);

   thrust::transform(d_x.begin(), d_x.end(), d_y.begin(), rotation_functor(R));

   thrust::host_vector<Vector> h_y = d_y;

   std::cout << "1: x: " << h_x[0] << ", y: " << h_y[0] << std::endl;
   std::cout << "2: x: " << h_x[1] << ", y: " << h_y[1] << std::endl;
   std::cout << "3: x: " << h_x[2] << ", y: " << h_y[2] << std::endl;
}
