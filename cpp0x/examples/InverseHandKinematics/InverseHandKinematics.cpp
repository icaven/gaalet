#include "gaalet.h"

#include <osgViewer/Viewer>
#include <osgGA/StateSetManipulator>
#include <osgGA/TrackballManipulator>
#include <osg/PositionAttitudeTransform>
#include <osg/ShapeDrawable>

typedef gaalet::algebra<gaalet::signature<3,0,1>> ma;

int main()
{
   ma::mv<0>::type one={1.0};
   ma::mv<1>::type e1={1.0};
   ma::mv<2>::type e2={1.0};
   ma::mv<4>::type e3={1.0};
   ma::mv<8>::type e0={1.0};
   ma::mv<0xf>::type I= (e1^e2^e3^e0);

   double K_d_r = 10.0;
   double K_d_t = 10.0;
   double K_1 = 1.0;
   double K_2 = 1.0;

   auto T_s = one + 0.5*I*(10.0*e2*e3 + 0.0*e3*e1 + -10.0*e1*e2);
   double phi_s = 0.0*M_PI;
   auto R_s = one*cos(-phi_s*0.5) + sin(-phi_s*0.5)*(0.0*e2*e3 + 1.0*e3*e1 + 0.0*e1*e2);
   auto D_s = eval(T_s*R_s);

   //Ball joint R_1:
   auto T_1 = eval(one + 0.5*I*(7.0*e2*e3 + 0.0*e3*e1 + 0.0*e1*e2));
   double phi_1 = -0.2*M_PI;
   auto R_1_0 = eval(one*cos(-phi_1*0.5) + sin(-phi_1*0.5)*(0.0*e2*e3 + 1.0*e3*e1 + 0.0*e1*e2));

   ma::mv<3,5,6>::type B_1;
   ma::mv<3,5,6>::type B_1_old;
   ma::mv<3,5,6>::type grad_B1_U_old;

   //Ball joint R_2:
   auto T_2 = eval(one + 0.5*I*(7.0*e2*e3 + 0.0*e3*e1 + 0.0*e1*e2));
   double phi_2 = 0.4*M_PI;
   auto R_2_0 = eval(one*cos(-phi_2*0.5) + sin(-phi_2*0.5)*(0.0*e2*e3 + 1.0*e3*e1 + 0.0*e1*e2));

   ma::mv<3,5,6>::type B_2;
   ma::mv<3,5,6>::type B_2_old;
   ma::mv<3,5,6>::type grad_B2_U_old;

   //Visualisation with OpenSceneGraph
   osg::Group* sceneRoot = new osg::Group;

   auto revconjD_s = D_s;
   revconjD_s[1] = -revconjD_s[1];
   revconjD_s[2] = -revconjD_s[2];
   revconjD_s[3] = -revconjD_s[3];
   auto X_s = eval(D_s*revconjD_s);

   osg::Sphere* sphere = new osg::Sphere(osg::Vec3(-X_s[4], -X_s[5], -X_s[6]), 0.3);
   osg::ShapeDrawable* sphereDrawable = new osg::ShapeDrawable(sphere);
   osg::Geode* sphereGeode = new osg::Geode();
   sphereGeode->addDrawable(sphereDrawable);
   sceneRoot->addChild(sphereGeode);

   osg::Box* rod1 = new osg::Box(osg::Vec3(3.5, 0.0, 0.0), 7.0, 0.2, 0.2);
   osg::ShapeDrawable* rod1Drawable = new osg::ShapeDrawable(rod1);
   osg::Geode* rod1Geode = new osg::Geode();
   rod1Geode->addDrawable(rod1Drawable);
   osg::PositionAttitudeTransform* rod1Transform = new osg::PositionAttitudeTransform();
   rod1Transform->addChild(rod1Geode);
   sceneRoot->addChild(rod1Transform);

   osg::Box* rod2 = new osg::Box(osg::Vec3(3.5, 0.0, 0.0), 7.0, 0.2, 0.2);
   osg::ShapeDrawable* rod2Drawable = new osg::ShapeDrawable(rod2);
   osg::Geode* rod2Geode = new osg::Geode();
   rod2Geode->addDrawable(rod2Drawable);
   osg::PositionAttitudeTransform* rod2Transform = new osg::PositionAttitudeTransform();
   rod2Transform->addChild(rod2Geode);
   sceneRoot->addChild(rod2Transform);

   osgViewer::Viewer viewer;
   viewer.setSceneData( sceneRoot );
   viewer.addEventHandler( new osgGA::StateSetManipulator(viewer.getCamera()->getOrCreateStateSet()) );
   if (!viewer.getCameraManipulator() && viewer.getCamera()->getAllowEventFocus())
   {
      viewer.setCameraManipulator(new osgGA::TrackballManipulator());
   }
   viewer.setReleaseContextAtEndOfFrameHint(false);

   if (!viewer.isRealized())
   {
      viewer.realize();
   }


   //double alpha_0 = 0.001;
   //for(int i=0; i<100; ++i) {
   double frameTime = 0.0;
   double sumFrameTime = 0.0;
   double minFrameTime = 0.1;
   double timer = 0.0;
   unsigned int counter = 0;
   while(!viewer.done()) {
      osg::Timer_t startFrameTick = osg::Timer::instance()->tick();


      auto R_1_d = (one+B_1)*!(one-B_1);
      auto R_2_d = (one+B_2)*!(one-B_2);
      double alpha = 0.001;

      auto D_d = eval(~T_1*~R_1_0*~R_1_d*~T_2*~R_2_0*~R_2_d*D_s);
      auto B_d_inv_part = eval(!(2.0*one+D_d+~D_d));
      auto B_d = eval((D_d-~D_d)*B_d_inv_part);

      //Ball joint R_1:
      auto grad_B1_Dd = eval(~T_1*~R_1_0*(-2.0*(one-B_1)*!(one+B_1)*!(one+B_1))*~T_2*~R_2_0*~R_2_d*D_s);
      auto grad_B1_revDd = eval(~D_s*R_2_d*R_2_0*T_2*(2.0*!(one-B_1))*R_1_0*T_1);
      auto grad_B1_Bd = (grad_B1_Dd-grad_B1_revDd)*B_d_inv_part-B_d*B_d_inv_part*(grad_B1_Dd+grad_B1_revDd);
      auto grad_B1_Ud = (grad_B1_Bd*(K_d_r*part<3,5,6>(B_d)+K_d_t*part<9,0xa,0xc>(B_d)));
     
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
      auto grad_B2_Dd = eval(~T_1*~R_1_0*~R_2_0*~T_2*(-2.0*(one-B_2)*!(one+B_2)*!(one+B_2))*~R_2_d*D_s);
      auto grad_B2_revDd = eval(~D_s*(2.0*!(one-B_2))*R_2_0*T_2*R_2_d*R_1_0*T_1);
      auto grad_B2_Bd = (grad_B2_Dd-grad_B2_revDd)*B_d_inv_part-B_d*B_d_inv_part*(grad_B2_Dd+grad_B2_revDd);
      auto grad_B2_Ud = (grad_B2_Bd*(K_d_r*part<3,5,6>(B_d)+K_d_t*part<9,0xa,0xc>(B_d)));
     
      auto grad_B2_U1 = K_2*B_2;

      auto grad_B2_U = grad_B2_U1 + grad_B2_Ud;

      auto delta_grad_B2_U = eval(grad_B2_U - grad_B2_U_old);
      grad_B2_U_old = grad_B2_U;
      auto delta_B2 = eval(B_2-B_2_old);
      B_2_old = B_2;

      B_2 = B_2 - alpha*grad_B2_U;


      auto D_1 = eval(R_1_d*R_1_0*T_1);
      auto D_2 = eval(R_2_d*R_2_0*T_2);
      auto D_12 = eval(D_1*D_2);
      auto revconjD_1 = D_1;
      revconjD_1[1] = -revconjD_1[1];
      revconjD_1[2] = -revconjD_1[2];
      revconjD_1[3] = -revconjD_1[3];
      auto X_1 = eval(D_1*revconjD_1);
      auto revconjD_12 = D_12;
      revconjD_12[1] = -revconjD_12[1];
      revconjD_12[2] = -revconjD_12[2];
      revconjD_12[3] = -revconjD_12[3];
      auto X_2 = eval(D_12*revconjD_12);
      //std::cout << "X_1: " << X_1 << ", X_2: " << X_2 << std::endl;
      rod1Transform->setAttitude(osg::Quat(-D_1[3], D_1[2], -D_1[1], D_1[0]));
      rod2Transform->setPosition(osg::Vec3d(-X_1[4], -X_1[5], -X_1[6]));
      rod2Transform->setAttitude(osg::Quat(-D_12[3], D_12[2], -D_12[1], D_12[0]));

      double E_Bd = 0.5*(K_d_r*B_d[0]*B_d[0]+K_d_r*B_d[1]*B_d[1]+K_d_r*B_d[2]*B_d[2]+K_d_t*B_d[3]*B_d[3]+K_d_t*B_d[4]*B_d[4]+K_d_t*B_d[5]*B_d[5]);
      std::cout << "i=" << counter << ": B_1: " << B_1 << ", B_2: " << B_2 << std::endl;
      std::cout << "\tB_d: " << B_d << ", E_Bd: " << E_Bd << std::endl;
      std::cout << "\tgrad_B1_U: " << grad_B1_U << std::endl;
      std::cout << "\tgrad_B2_U: " << grad_B2_U << std::endl;

      viewer.frame();

      //work out if we need to force a sleep to hold back the frame rate
      osg::Timer_t endFrameTick = osg::Timer::instance()->tick();
      frameTime = osg::Timer::instance()->delta_s(startFrameTick, endFrameTick);

      sumFrameTime += frameTime;
      if(counter==1000) {
         std::cout << "Average frame time: " << sumFrameTime/1000.0 << std::endl;
         sumFrameTime = 0.0;
         counter = 0;
      }
      else {
         counter++;
      }

      timer += frameTime;

      if (frameTime < minFrameTime) OpenThreads::Thread::microSleep(static_cast<unsigned int>(1000000.0*(minFrameTime-frameTime)));

   }

}
