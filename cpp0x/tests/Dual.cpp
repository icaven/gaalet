#include "gaalet.h"

typedef gaalet::algebra< gaalet::signature<3,0> > em;

typedef gaalet::algebra< gaalet::signature<0,6,2> > sa;

typedef gaalet::algebra< gaalet::signature<3, 0, 1> > pga;


int main()
{
   em::mv<7>::type Ie = {1.0};
   em::mv<1,2,4>::type a = {1.0,2.0,3.0};

   std::cout << "a: " << a << ", *a: " << dual(a) << ", **a: " << dual(dual(a)) << std::endl;
   std::cout << "I: " << Ie << ", *I: " << dual(Ie) << ", **I: " << dual(dual(Ie)) << std::endl;
   
   std::cout << "a*~I: " << a*(~Ie) << ", a*~I*I: " << a*(~Ie)*Ie << std::endl;
   std::cout << "a*~I: " << a*(~Ie) << ", a*~I*~I: " << a*(~Ie)*(~Ie) << std::endl;

   static const sa::mv<1>::type e1={1.0};
   static const sa::mv<2>::type e2={1.0};
   static const sa::mv<4>::type e3={1.0};
   static const sa::mv<0x40>::type e0={1.0};
   static const sa::mv<(1<<(6+2))-1>::type I={1.0};  // The pseudoscalar for the sa algebra

   std::cout << "I: " << I << ", *I: " << dual(I) << ", **I: " << dual(dual(I)) << std::endl;

   std::cout << "dual(e1*e2): " << dual(e1*e2) << std::endl;
   std::cout << "dual(e3*e0): " << dual(e3*e0) << std::endl;
   std::cout << "*(*(e1*e2) ^ *(e3*e0)): " << dual(dual(e1*e2) ^ dual(e3*e0)) << std::endl;
   
   static const pga::mv<0x0>::type one = { 1.0 };
   static const pga::mv<1>::type e1_pga={1.0};
   static const pga::mv<2>::type e2_pga={1.0};
   static const pga::mv<4>::type e3_pga={1.0};
   static const pga::mv<8>::type e0_pga={1.0};
   static const pga::mv<(1<<(3+1))-1>::type Ipga={1.0};  // The pseudoscalar for the pga algebra

   std::cout << "one: " << one << ", *one: " << dual(one) << ", **Ipga: " << dual(dual(one)) << std::endl;
   std::cout << "one: " << one << ", *one: " << ::dual(one) << ", **Ipga: " << ::dual(::dual(one)) << std::endl;
   std::cout << "Ipga: " << Ipga << ", *Ipga: " << dual(Ipga) << ", **Ipga: " << dual(dual(Ipga)) << std::endl;
   std::cout << "e0_pga^e1_pga: " << (e0_pga^e1_pga) << std::endl;
   std::cout << "e2_pga^e3_pga: " << (e2_pga^e3_pga) << std::endl;

   std::cout << "dual(e0_pga^e1_pga): " << dual(e0_pga^e1_pga) << std::endl;
   std::cout << "dual(e2_pga^e3_pga): " << dual(e2_pga^e3_pga) << std::endl;

}
