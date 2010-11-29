#ifndef __CAR_DYNAMICS_H
#define __CAR_DYNAMICS_H

#include "gaalet.h"
#include "MagicFormula2004.h"
#include <tuple>

namespace cardyn {


//defintion of basisvectors, null basis, pseudoscalars, helper unit scalar
typedef gaalet::algebra< gaalet::signature<4,1> > cm;
cm::mv<0x01>::type e1 = {1.0};
cm::mv<0x02>::type e2 = {1.0};
cm::mv<0x04>::type e3 = {1.0};
cm::mv<0x08>::type ep = {1.0};
cm::mv<0x10>::type em = {1.0};

cm::mv<0x00>::type one = {1.0};

cm::mv<0x08, 0x10>::type e0 = 0.5*(em-ep);
cm::mv<0x08, 0x10>::type einf = em+ep;

cm::mv<0x18>::type E = ep*em;

cm::mv<0x1f>::type Ic = e1*e2*e3*ep*em;
cm::mv<0x07>::type Ie = e1*e2*e3;

typedef cm::mv<0x03, 0x05, 0x06, 0x09, 0x0a, 0x0c, 0x11, 0x12, 0x14>::type S_type;
typedef cm::mv<0x00, 0x03, 0x05, 0x06, 0x09, 0x0a, 0x0c, 0x0f, 0x11, 0x12, 0x14, 0x17>::type D_type;


typedef std::tuple<
         double,    //Vertical wheel distance to ground, front lef
         double,    //Vertical wheel distance to ground, front right
         double,    //Vertical wheel distance to ground, rear left
         double,    //Vertical wheel distance to ground, rear right
         double,    //Steering angle
         double,    //Transmission
         double,    //Gas pedal position [0.0-1.0]
         double,    //Brake shoe force
         double     //Clutch coefficient
        > InputVector;

typedef std::tuple<
      D_type,     //Body displacement
      S_type,     //Screw velocity
      double,   
      double,
      double,                  //Spring damper compression, wheel front left
      double,
      double,                  //Spring damper compression, wheel front right
      double,
      double,                  //Spring damper compression, wheel rear left
      double,
      double,                  //Spring damper compression, wheel rear right
      double,
      double,                  //Wheel angular velocity, wheel front left
      double,                  //Wheel angular velocity, wheel front right
      double,                  //Wheel angular velocity, wheel rear left
      double,                  //Wheel angular velocity, wheel rear right
      double                   //Engine speed
   > StateVector;

struct StateEquation
{
   StateEquation( const InputVector& input_,
                  const magicformula2004::ContactWrench& tyre_fl_,
                  const magicformula2004::ContactWrench& tyre_fr_,
                  const magicformula2004::ContactWrench& tyre_rl_,
                  const magicformula2004::ContactWrench& tyre_rr_)
      :  input(input_),
         tyre_fl(tyre_fl_),
         tyre_fr(tyre_fr_),
         tyre_rl(tyre_rl_),
         tyre_rr(tyre_rr_)
   {
      g[0] = -9.81;

      cm::mv<1,2,4>::type x_wfl = {1.410, 0.747, -0.4};
      r_wfl = x_wfl + 0.5*(x_wfl&x_wfl)*einf + e0;
      cm::mv<1,2,4>::type x_wfr = {1.410, -0.747, -0.4};
      r_wfr = x_wfr + 0.5*(x_wfr&x_wfr)*einf + e0;
      cm::mv<1,2,4>::type x_wrl = {-0.940, 0.812, -0.4};
      r_wrl = x_wrl + 0.5*(x_wrl&x_wrl)*einf + e0;
      cm::mv<1,2,4>::type x_wrr = {-0.940, -0.812, -0.4};
      r_wrr = x_wrr + 0.5*(x_wrr&x_wrr)*einf + e0;

      q_w[0] = cos(0.5*M_PI); q_w[1] = sin(0.5*M_PI);
     
      i_g.resize(7);
      i_g[0] = -3.6;
	   i_g[1] = 0.0;
	   i_g[2] = 3.6;
	   i_g[3] = 2.19;
   	i_g[4] = 1.41;
	   i_g[5] = 1.0;
	   i_g[6] = 0.83; 

      R_n_wfl = exp((-0.5)*Ie*(M_PI*0.05*e1 + M_PI*0.1*e2 + 0.0*e3));
      R_n_wfr = exp((-0.5)*Ie*(-M_PI*0.05*e1 + M_PI*0.1*e2 + 0.0*e3));
      R_n_wrl = exp((-0.5)*Ie*(M_PI*0.05*e1 + 0.0*e2 + 0.0*e3));
      R_n_wrr = exp((-0.5)*Ie*(-M_PI*0.05*e1 + 0.0*e2 + 0.0*e3));
   }

   StateVector operator()(const double& t, const StateVector& oldState) const
   {
      //state
      const auto& D_b = std::get<0>(oldState);
      //const auto& p_b = std::get<0>(oldState);
      //const auto& dp_b = std::get<1>(oldState);
      const auto& V_b = std::get<1>(oldState);
      //const auto& q_b = std::get<2>(oldState);
      //const auto& w_b = std::get<3>(oldState);
      const auto& u_wfl= std::get<4>(oldState);
      const auto& du_wfl= std::get<5>(oldState);
      const auto& u_wfr= std::get<6>(oldState);
      const auto& du_wfr= std::get<7>(oldState);
      const auto& u_wrl= std::get<8>(oldState);
      const auto& du_wrl= std::get<9>(oldState);
      const auto& u_wrr= std::get<10>(oldState);
      const auto& du_wrr= std::get<11>(oldState);
      const auto& w_wfl= std::get<12>(oldState);
      const auto& w_wfr= std::get<13>(oldState);
      const auto& w_wrl= std::get<14>(oldState);
      const auto& w_wrr= std::get<15>(oldState);
      const auto& w_e= std::get<16>(oldState);

      //input
      const double& Dv_wfl = std::get<0>(input);
      const double& Dv_wfr = std::get<1>(input);
      const double& Dv_wrl = std::get<2>(input);
      const double& Dv_wrr = std::get<3>(input);
      const double& steerAngle = std::get<4>(input);
      const double& i_pt = std::get<5>(input);
      const double& s_gp = std::get<6>(input);
      const double& F_b = std::get<7>(input);


      //Ackermann steering
      double cotSteerAngle = (r_wfl[0]-r_wrl[0])*(1.0/tan(steerAngle));
      double angleFL = atan(1.0/(cotSteerAngle - w_wn/(v_wn*2.0)));
      double angleFR = atan(1.0/(cotSteerAngle + w_wn/(v_wn*2.0)));
      gaalet::mv<0,3>::type q_wfl = {cos(angleFL*0.5), sin(angleFL*0.5)};
      gaalet::mv<0,3>::type q_wfr = {cos(angleFR*0.5), sin(angleFR*0.5)};

      //wheel velocity in body frame:
      auto dr_wfl = grade<1>(q_wfl*((V_b&r_wfl)-du_wfl*e3)*(!q_wfl));
      auto dr_wfr = grade<1>(q_wfr*((V_b&r_wfr)-du_wfr*e3)*(!q_wfr));
      auto dr_wrl = grade<1>((V_b&r_wrl)-du_wrl*e3);
      auto dr_wrr = grade<1>((V_b&r_wrr)-du_wrr*e3);

      //wheel rotors:
      auto R_wfl = R_n_wfl*exp(e2*e3*u_wfl*(-0.5));
      auto R_wfr = R_n_wfr*exp(e2*e3*u_wfr*(0.5));
      auto R_wrl = R_n_wrl*exp(e2*e3*u_wrl*(-0.5));
      auto R_wrr = R_n_wrr*exp(e2*e3*u_wrr*(0.5));
      
      //Suspension spring damper force:
      auto Fsd_wfl = (u_wfl*k_wf+du_wfl*d_wf)*(-1.0);
      auto Fsd_wfr = (u_wfr*k_wf+du_wfr*d_wf)*(-1.0);
      auto Fsd_wrl = (u_wrl*k_wr+du_wrl*d_wr)*(-1.0);
      auto Fsd_wrr = (u_wrr*k_wr+du_wrr*d_wr)*(-1.0);

      //Tyre forces and moments:
      auto W_wfl = ((!q_w)*tyre_fl( Dv_wfl, //distance difference with respect to camber angle?
            R_wfl,
            part<1,2,4,5>(q_w*(dr_wfl + w_wfl*e1*e3)*(!q_w)))*q_w);
      auto W_wfr = ((!q_w)*tyre_fr( Dv_wfr, //distance difference with respect to camber angle?
            R_wfr,
            part<1,2,4,5>(q_w*(dr_wfr + w_wfr*e1*e3)*(!q_w)))*q_w);
      auto W_wrl = ((!q_w)*tyre_rl( Dv_wrl, //distance difference with respect to camber angle?
            R_wrl,
            part<1,2,4,5>(q_w*(dr_wrl + w_wrl*e1*e3)*(!q_w)))*q_w);
      auto W_wrr = ((!q_w)*tyre_rr( Dv_wrr, //distance difference with respect to camber angle?
            R_wrr,
            part<1,2,4,5>(q_w*(dr_wrr + w_wrr*e1*e3)*(!q_w)))*q_w);

      //Body acceleration:
      auto ddp_b_b = eval(grade<1>((((grade<1>((!q_wfl)*part<1, 2>(W_wfl)*q_wfl+(!q_wfr)*part<1, 2>(W_wfr)*q_wfr+part<1, 2>(W_wrl)+part<1, 2>(W_wrr))+(Fsd_wfl+Fsd_wfr+Fsd_wrl+Fsd_wrr)*e3)*(1.0/m_b)))) + grade<1>((!part<0,3,5,6>(D_b))*g*part<0,3,5,6>(D_b)));

      auto w_b_b = eval(Ie&V_b);
      double k_arb = this->k_arb;
      cm::mv<1,2,4>::type t_b_b = (-1.0)*Ie*((part<1,2,4>(r_wfl)^(Fsd_wfl*e3+grade<1>((!q_wfl)*part<1,2>(W_wfl)*q_wfl)-(u_wfl-u_wfr)*e3*k_arb)) + (part<1,2,4>(r_wfr)^(Fsd_wfr*e3+grade<1>((!q_wfr)*part<1,2>(W_wfr)*q_wfr)+(u_wfl-u_wfr)*e3*k_arb)) + (part<1,2,4>(r_wrl)^(Fsd_wrl*e3+part<1,2>(W_wrl))) + (part<1,2,4>(r_wrr)^(Fsd_wrr*e3+part<1,2>(W_wrr))));
      //cm::mv<1,2,4>::type t_b_b = (-1.0)*Ie*((part<1,2,4>(r_wfl)^(Fsd_wfl*z-(u_wfl-u_wfr)*z*k_arb)) + (part<1,2,4>(r_wfr)^(Fsd_wfr*z+(u_wfl-u_wfr)*z*k_arb)) + (part<1,2,4>(r_wrl)^(Fsd_wrl*z)) + (part<1,2,4>(r_wrr)^(Fsd_wrr*z+part<1,2>(W_wrr))));
      //cm::mv<1,2,4>::type t_b_b = (-1.0)*Ie*((part<1,2,4>(r_wfl)^(Fsd_wfl*z)) + (part<1,2,4>(r_wfr)^(Fsd_wfr*z)) + (part<1,2,4>(r_wrl)^(Fsd_wrl*z)) + (part<1,2,4>(r_wrr)^(Fsd_wrr*z)));
      cm::mv<1,2,4>::type dw_b_b;
      double In_1 = 590.0, In_2 = 1730.0, In_3 = 1950.0;
      dw_b_b[0] = (t_b_b[0] - (In_3-In_2)*w_b_b[1]*w_b_b[2])/In_1;
      dw_b_b[1] = (t_b_b[1] - (In_1-In_3)*w_b_b[2]*w_b_b[0])/In_2;
      dw_b_b[2] = (t_b_b[2] - (In_2-In_1)*w_b_b[0]*w_b_b[1])/In_3;

      auto dV_b = (-1.0)*Ie*dw_b_b + einf*ddp_b_b;

      StateVector newState(
         part_type<D_type>(D_b*V_b*0.5),
         dV_b,
         0.0,
         0.0,
         du_wfl,
         (Fsd_wfl - W_wfl.element<4>())*(1.0/m_w),
         du_wfr,
         (Fsd_wfr - W_wfr.element<4>())*(1.0/m_w),
         du_wrl,
         (Fsd_wrl - W_wrl.element<4>())*(1.0/m_w),
         du_wrr,
         (Fsd_wrr - W_wrr.element<4>())*(1.0/m_w),
         (W_wfl.element<1>()*(-r_w) - tanh(w_wfl*d_b)*mu_b*F_b)*(1.0/I_w),
         (W_wfr.element<1>()*(-r_w) - tanh(w_wfr*d_b)*mu_b*F_b)*(1.0/I_w),
         (W_wrl.element<1>()*(-r_w) - tanh(w_wrl*d_b)*mu_b*F_b + i_pt*(w_e-(w_wrl+w_wrr)*i_pt*0.5))*(1.0/I_w),
         (W_wrr.element<1>()*(-r_w) - tanh(w_wrr*d_b)*mu_b*F_b + i_pt*(w_e-(w_wrl+w_wrr)*i_pt*0.5))*(1.0/I_w),
         (s_gp*(w_e*w_e*a_e + w_e*b_e + c_e) - (w_e-(w_wrl+w_wrr)*i_pt*0.5) - w_e*d_e)*(1.0/I_e)
      );

      return std::move(newState);
   }


   const InputVector& input;
   magicformula2004::ContactWrench tyre_fl;
   magicformula2004::ContactWrench tyre_fr;
   magicformula2004::ContactWrench tyre_rl;
   magicformula2004::ContactWrench tyre_rr;

   cm::mv<4>::type g;

   //Wheel positions in car body frame
   cm::mv<1,2,4,8,0x10>::type r_wfl;
   cm::mv<1,2,4,8,0x10>::type r_wfr;
   cm::mv<1,2,4,8,0x10>::type r_wrl;
   cm::mv<1,2,4,8,0x10>::type r_wrr;

   cm::mv<0,6>::type q_w;

   //Carbody
   static const double m_b = 1450.0;
   static const double r_b = 3.0;
   //Sphere
   //static const double I_b = 2.0/5.0*m_b*r_b*r_b;

   //Wheel
   static const double m_w = 20;
   static const double r_w = 0.325;
   static const double I_w = 2.3;
   static const double u_wn = 0.4;
   static const double v_wn = 1.3;
   static const double w_wn = 0.7;
   static const double k_wf = 17400.0;
   static const double k_wr = 26100.0;
   static const double d_wf = 2600.0;
   static const double d_wr = 2600.0;
   gaalet::mv<0,3,5,6>::type R_n_wfl;
   gaalet::mv<0,3,5,6>::type R_n_wfr;
   gaalet::mv<0,3,5,6>::type R_n_wrl;
   gaalet::mv<0,3,5,6>::type R_n_wrr;

   //Braking system
   static const double mu_b = 0.135;
   static const double d_b = 0.01;

   //Anti roll bar
   //static const double k_arb = 50000;
   static const double k_arb = 100000;

   //Clutch
   static const double k_cn = 1.5;

   //Engine
   static const double a_e = -0.000862;
   static const double b_e = 0.83;
   static const double c_e = 400;
   static const double I_e = 0.5;
   static const double d_e = 0.5;

   //Transmission
   std::vector<double> i_g;
   static const double i_a = 3.5;
};

}  //end namespace cardyn

#endif
