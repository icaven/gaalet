#include "gaalet.h"
#include <iostream>
#include <sys/time.h>

template <class P, class EInf>
auto sphere(P& center, double radius, EInf& einf) {
    return dual(center - radius*radius * 0.5 * einf);
}

template <class P, class EInf>
auto plane(P& v, double h, EInf& einf) {
    return dual(v - h * einf);
}

template <class A, class B>
auto vee(A& a, B& b) {
    return ::dual(::dual(a) ^ ::dual(b));
}

template <class P, class EO, class EInf>
auto up(P& p, EO& e0, EInf& einf) {
    return p + e0 + 0.5*p*p*einf;
}

int main()
{
   gaalet::cm::mv<0x01>::type e1 = {1.0};
   gaalet::cm::mv<0x02>::type e2 = {1.0};
   gaalet::cm::mv<0x04>::type e3 = {1.0};
   gaalet::cm::mv<0x08>::type ep = {1.0};
   gaalet::cm::mv<0x10>::type em = {1.0};

   gaalet::cm::mv<0x00>::type one = {1.0};

   gaalet::cm::mv<0x08, 0x10>::type e0 = 0.5*(em-ep);
   //auto e0 = 0.5*(em-ep);
   gaalet::cm::mv<0x08, 0x10>::type einf = em+ep;
   //auto einf = em+ep;

   double x = 2.0;
   double y = 1.0;
   double z = 0.5;

   auto a_point = e1*x + e2*y + e3*z;
   auto a_conformal_point = up(a_point, e0, einf);
   auto P = eval(e1*x + e2*y + e3*z + e0*1.0 + einf*(x*x+y*y+z*z)*0.5);

   gaalet::cm::mv<0x08, 0x10>::type S;
   gaalet::mv<0x09, 0x0a, 0x0c, 0x11, 0x12, 0x14, 0x18>::type C;

   timeval start, end;

   double r = 0.5;
   
   
   auto cS = sphere(e0, r, einf);
   auto cP = plane(a_point, -0.5*eval(magnitude2(a_point)).template element<0>(), einf);
   
   std::cout << "Squared magnitude of a_point " << magnitude2(a_point) << std::endl;
   std::cout << "(x*x+y*y+z*z) " << (x*x+y*y+z*z) << std::endl;

   gettimeofday(&start, 0);
   for(int i=0; i<100000000; ++i, r+=0.000001) {
      S = e0 - einf*0.5*r*r;
      C = (S^(P+(P&S)*einf));
   }
   gettimeofday(&end, 0);
   std::cout << "Time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)*1e-6 << ", circle C: " << C << std::endl;

   gettimeofday(&start, 0);
   for(int i=0; i<100000000; ++i, r+=0.000001) {
      auto cS = sphere(e0, r, einf);
      C = vee(cS, cP);
   }
   gettimeofday(&end, 0);
   std::cout << "Time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)*1e-6 << ", circle2 C: " << C << std::endl;

   double C_opt[17];

   r = 0.5;
   gettimeofday(&start, 0);
   for(int i=0; i<100000000; ++i, r+=0.000001) {
      C_opt[9]=0.5*x*r*r;
      C_opt[10]=-1*x;
      C_opt[12]=0.5*y*r*r;
      C_opt[13]=-1*y;
      C_opt[14]=0.5*z*r*r;
      C_opt[15]=-1*z;
      C_opt[16]=-1*r*r;
   }
   gettimeofday(&end, 0);

   std::cout << "Time: " << (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)*1e-6 
             << ", C: [" << C_opt[9] << ", " << C_opt[10] << ", " 
             << C_opt[12] << ", " << C_opt[13] << ", " << C_opt[14] << ", " 
             << C_opt[15] << ", " << C_opt[16] << "]" << std::endl;
}
