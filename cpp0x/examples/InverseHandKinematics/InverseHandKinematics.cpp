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

   auto T_s = one + 0.5*I*(10.0*e2*e3 + 0.0*e3*e1 + 0.0*e1*e2);
   double phi_s = 0.25*M_PI;
   auto R_s = one*cos(phi_s*0.5) + sin(phi_s*0.5)*(0.0*e2*e3 + 1.0*e3*e1 + 0.0*e1*e2);
   auto D_s = eval(T_s*R_s);

   //Ball joint R_1:
   auto T_1 = eval(one + 0.5*I*(7.0*e2*e3 + 0.0*e3*e1 + 0.0*e1*e2));
   double phi_1 = -0.25*M_PI;
   auto R_1_0 = eval(one*cos(phi_1*0.5) + sin(phi_1*0.5)*(0.0*e2*e3 + 1.0*e3*e1 + 0.0*e1*e2));

   ma::mv<3,5,6>::type B_1;
   ma::mv<3,5,6>::type B_1_old;
   ma::mv<3,5,6>::type grad_B1_U_old;

   //Ball joint R_2:
   auto T_2 = eval(one + 0.5*I*(7.0*e2*e3 + 0.0*e3*e1 + 0.0*e1*e2));
   double phi_2 = 0.25*M_PI;
   auto R_2_0 = eval(one*cos(phi_2*0.5) + sin(phi_2*0.5)*(0.0*e2*e3 + 1.0*e3*e1 + 0.0*e1*e2));

   ma::mv<3,5,6>::type B_2;
   ma::mv<3,5,6>::type B_2_old;
   ma::mv<3,5,6>::type grad_B2_U_old;

   //double alpha_0 = 0.001;
   for(int i=0; i<100; ++i) {
      auto R_1_d = (one+B_1)*!(one-B_1);
      auto R_2_d = (one+B_2)*!(one-B_2);
      std::cout << "i=" << i << ": B_1: " << B_1 << ", B_2: " << B_2 << std::endl;
      double alpha = 0.001;

      auto D_d = eval(~T_1*~R_1_0*~R_1_d*~T_2*~R_2_0*~R_2_d*D_s);
      auto B_d = (D_d-~D_d)*~(2.0*one+D_d+~D_d);

      //Ball joint R_1:
      auto grad_B1_Dd = ~T_1*~R_1_0*(-2.0*(one-B_1)*!(one+B_1)*!(one+B_1))*~T_2*~R_2_0*~R_2_d*D_s;
      auto grad_B1_revDd = ~D_s*R_2_d*R_2_0*T_2*(2.0*!(one-B_1))*R_1_0*T_1;
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

      B_1 = B_1 - alpha*grad_B1_U;

      //Ball joint R_2:
      auto grad_B2_Dd = ~T_1*~R_1_0*~R_2_0*~T_2*(-2.0*(one-B_2)*!(one+B_2)*!(one+B_2))*~R_2_d*D_s;
      auto grad_B2_revDd = ~D_s*(2.0*!(one-B_2))*R_2_0*T_2*R_2_d*R_1_0*T_1;
      auto grad_B2_Bd = (grad_B2_Dd-grad_B2_revDd)*!(2.0*one+D_d+~D_d)-(D_d-~D_d)*!(2.0*one+D_d+~D_d)*!(2.0*one+D_d+~D_d)*(grad_B2_Dd+grad_B2_revDd);
      auto grad_B2_Ud = (grad_B2_Bd*(K_d*B_d));
     
      auto grad_B2_U1 = K_1*B_2;

      auto grad_B2_U = grad_B2_U1 + grad_B2_Ud;

      auto delta_grad_B2_U = eval(grad_B2_U - grad_B2_U_old);
      grad_B2_U_old = grad_B2_U;
      auto delta_B2 = eval(B_2-B_2_old);
      B_2_old = B_2;

      B_2 = B_2 - alpha*grad_B2_U;

   }

   auto R_1_d = (one+B_1)*!(one-B_1);
   auto R_2_d = (one+B_2)*!(one-B_2);
   std::cout << "R_s: " << R_s << ", R_1_d*R_1_0*R_2_d*R_2_0: " << R_1_d*R_1_0*R_2_d*R_2_0 << std::endl;
}
