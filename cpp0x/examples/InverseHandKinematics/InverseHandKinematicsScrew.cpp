#include "gaalet.h"

#include <osgViewer/Viewer>
#include <osgGA/StateSetManipulator>
#include <osgGA/TrackballManipulator>
#include <osg/PositionAttitudeTransform>
#include <osg/ShapeDrawable>

typedef gaalet::algebra<gaalet::signature<3,0,1>> ma;
typedef ma::mv<0, 3, 5, 6, 9, 0xa, 0xc>::type Motor;

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

   auto T_t = one + 0.5*I*(10.0*e2*e3 + 0.0*e3*e1 + -10.0*e1*e2);
   double phi_s = 0.0*M_PI;
   auto R_t = one*cos(-phi_s*0.5) + sin(-phi_s*0.5)*(0.0*e2*e3 + 1.0*e3*e1 + 0.0*e1*e2);
   auto M_t = eval(T_t*R_t);

   //Ball joint R_1:
   auto T_1 = eval(one + 0.5*I*(7.0*e2*e3 + 0.0*e3*e1 + 0.0*e1*e2));
   double phi_1 = -0.2*M_PI;
   auto R_1_0 = eval(one*cos(-phi_1*0.5) + sin(-phi_1*0.5)*(0.0*e2*e3 + 1.0*e3*e1 + 0.0*e1*e2));

   ma::mv<3,5,6>::type S_1;

   //Ball joint R_2:
   auto T_2 = eval(one + 0.5*I*(7.0*e2*e3 + 0.0*e3*e1 + 0.0*e1*e2));
   double phi_2 = 0.4*M_PI;
   auto R_2_0 = eval(one*cos(-phi_2*0.5) + sin(-phi_2*0.5)*(0.0*e2*e3 + 1.0*e3*e1 + 0.0*e1*e2));

   ma::mv<3,5,6>::type S_2;

   //Visualisation with OpenSceneGraph
   osg::Group* sceneRoot = new osg::Group;

   auto revconjM_t = M_t;
   revconjM_t[1] = -revconjM_t[1];
   revconjM_t[2] = -revconjM_t[2];
   revconjM_t[3] = -revconjM_t[3];
   auto X_s = eval(M_t*revconjM_t);

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


      auto R_1_d = exp(0.5*S_1);
      auto R_2_d = exp(0.5*S_2);
      double alpha = 0.001;

      auto M_d = part_type<Motor>(~T_2*~R_2_0*~R_2_d*~T_1*~R_1_0*~R_1_d*M_t);
      auto S_d = eval(2.0*log(M_d));

      //Ball joint R_1:
      auto grad_S1_U1 = K_1*S_1;
      auto grad_S1_Ud = -1.0*(K_d_r*part<3,5,6>(S_d)+K_d_t*part<9,0xa,0xc>(S_d));
      auto grad_S1_U = grad_S1_U1 + grad_S1_Ud;

      S_1 = S_1 - alpha*grad_S1_U;

      //Ball joint R_2:
      auto grad_S2_U1 = K_1*S_2;
      auto grad_S2_Ud = -1.0*(K_d_r*part<3,5,6>(S_d)+K_d_t*part<9,0xa,0xc>(S_d));
      auto grad_S2_U = grad_S2_U1 + grad_S2_Ud;

      S_2 = S_2 - alpha*grad_S2_U;


      auto M_1 = eval(R_1_d*R_1_0*T_1);
      auto M_2 = eval(R_2_d*R_2_0*T_2);
      auto M_12 = eval(M_1*M_2);
      auto revconjM_1 = M_1;
      revconjM_1[1] = -revconjM_1[1];
      revconjM_1[2] = -revconjM_1[2];
      revconjM_1[3] = -revconjM_1[3];
      auto X_1 = eval(M_1*revconjM_1);
      auto revconjM_12 = M_12;
      revconjM_12[1] = -revconjM_12[1];
      revconjM_12[2] = -revconjM_12[2];
      revconjM_12[3] = -revconjM_12[3];
      auto X_2 = eval(M_12*revconjM_12);
      //std::cout << "X_1: " << X_1 << ", X_2: " << X_2 << std::endl;
      rod1Transform->setAttitude(osg::Quat(-M_1[3], M_1[2], -M_1[1], M_1[0]));
      rod2Transform->setPosition(osg::Vec3d(-X_1[4], -X_1[5], -X_1[6]));
      rod2Transform->setAttitude(osg::Quat(-M_12[3], M_12[2], -M_12[1], M_12[0]));

      double U_Sd = 0.5*(K_d_r*S_d[0]*S_d[0]+K_d_r*S_d[1]*S_d[1]+K_d_r*S_d[2]*S_d[2]+K_d_t*S_d[3]*S_d[3]+K_d_t*S_d[4]*S_d[4]+K_d_t*S_d[5]*S_d[5]);
      std::cout << "i=" << counter << ": S_1: " << S_1 << ", S_2: " << S_2 << std::endl;
      std::cout << "\tS_d: " << S_d << ", U_Sd: " << U_Sd << std::endl;
      std::cout << "\tgrad_S1_U: " << grad_S1_U << std::endl;
      std::cout << "\tgrad_S2_U: " << grad_S2_U << std::endl;

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
