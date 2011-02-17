#include "gaalet.h"

typedef gaalet::algebra<gaalet::signature<3,0,1>> ma;

int main()
{
   ma::mv<0>::type one={1.0};
   ma::mv<1>::type e1={1.0};
   ma::mv<2>::type e2={1.0};
   ma::mv<4>::type e3={1.0};
   ma::mv<8>::type e0={1.0};
   ma::mv<0xf>::type I= (e1^e2^e3^e0);

   double K_d = 10.0;
   double K_1 = 1.0;

   auto D_s = eval(one + 0.5*I*(10.0*e2*e3 + 0.0*e3*e1 + 0.0*e1*e2));

   auto T_1 = eval(one + 0.5*I*(10.0*e2*e3 + 0.0*e3*e1 + 0.0*e1*e2));
   double phi = 0.25*M_PI;
   auto R_1_0 = eval(one*cos(phi*0.5) + sin(phi*0.5)*(0.0*e2*e3 + 1.0*e3*e1 + 0.0*e1*e2));

   ma::mv<3,5,6>::type B_1;
   ma::mv<3,5,6>::type B_1_old;
   ma::mv<3,5,6>::type grad_B1_U_old;

   //double alpha_0 = 0.001;
   for(int i=0; i<100; ++i) {
      auto R_1_d = (one+B_1)*!(one-B_1);
      std::cout << "i=" << i << ": R_1_d: " << R_1_d << ", B_1: " << B_1 << std::endl;

      auto D_d = ~T_1*~R_1_0*~R_1_d*D_s;
      auto B_d = (D_d-~D_d)*~(2.0*one+D_d+~D_d);

      auto grad_B1_Dd = ~T_1*~R_1_0*(-2.0*(one-B_1)*!(one+B_1)*!(one+B_1))*D_s;
      auto grad_B1_revDd = ~D_s*(2.0*!(one-B_1))*R_1_0*T_1;
      auto grad_B1_Bd = (grad_B1_Dd-grad_B1_revDd)*!(2.0*one+D_d+~D_d)-(D_d-~D_d)*!(2.0*one+D_d+~D_d)*!(2.0*one+D_d+~D_d)*(grad_B1_Dd+grad_B1_revDd);
      auto grad_B1_Ud = (grad_B1_Bd*(K_d*B_d));
     
      auto grad_B1_U1 = K_1*B_1;

      auto grad_B1_U = grad_B1_U1 + grad_B1_Ud;

      auto delta_grad_B1_U = eval(grad_B1_U - grad_B1_U_old);
      grad_B1_U_old = grad_B1_U;
      auto delta_B1 = eval(B_1-B_1_old);
      B_1_old = B_1;
      //double alpha_two_step = eval(grade<0>(delta_B1&delta_grad_B1_U) * ~grade<0>(delta_grad_B1_U&delta_grad_B1_U));
      //double alpha_two_step = eval(grade<0>(delta_B1&delta_B1) * ~grade<0>(delta_B1&delta_grad_B1_U));
      //std::cout << "alpha: " << alpha_two_step << std::endl;
      double alpha = 0.001;

      B_1 = B_1 - alpha*grad_B1_U;
   }

   auto R_1_d = (one+B_1)*!(one-B_1);
   std::cout << "R_1_d: " << R_1_d << ", R_1_0: " << R_1_0 << ", R_1_d*R_1_0: " << R_1_d*R_1_0 << std::endl;
}
