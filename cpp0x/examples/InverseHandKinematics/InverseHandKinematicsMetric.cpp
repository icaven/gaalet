#include "gaalet.h"

#include <osgViewer/Viewer>
#include <osgGA/StateSetManipulator>
#include <osgGA/TrackballManipulator>
#include <osg/PositionAttitudeTransform>
#include <osg/ShapeDrawable>

typedef gaalet::algebra<gaalet::signature<3,0>> em;
typedef em::mv<1,2,4>::type Vector;
typedef em::mv<3,5,6>::type Bivector;
typedef em::mv<0, 3, 5, 6>::type Rotor;

int main()
{
   em::mv<0>::type one={1.0};
   em::mv<1>::type e1={1.0};
   em::mv<2>::type e2={1.0};
   em::mv<4>::type e3={1.0};
   auto I = eval(e1^e2^e3);

   double K_1 = 1.0;
   double K_2 = 1.0;

   auto t_t = eval(7.0*e1 + 0.0*e2 + -7.0*e3);
   double phi_t = 0.0*M_PI;
   auto R_t = one*cos(-phi_t*0.5) + sin(-phi_t*0.5)*(0.0*e2*e3 + 1.0*e3*e1 + 0.0*e1*e2);

   //Ball joint R_1:
   auto t_1 = eval(7.0*e1 + 0.0*e2 + 0.0*e3);
   double phi_1 = -0.2*M_PI;
   auto R_1_0 = eval(one*cos(-phi_1*0.5) + sin(-phi_1*0.5)*(0.0*e2*e3 + 1.0*e3*e1 + 0.0*e1*e2));

   Bivector m_1;

   //Ball joint R_2:
   auto t_2 = eval(7.0*e1 + 0.0*e2 + 0.0*e3);
   double phi_2 = 0.4*M_PI;
   auto R_2_0 = eval(one*cos(-phi_2*0.5) + sin(-phi_2*0.5)*(0.0*e2*e3 + 1.0*e3*e1 + 0.0*e1*e2));

   Bivector m_2;

   double lambda_t = 0.0;
   Bivector lambda_R;

   //Visualisation with OpenSceneGraph
   osg::Group* sceneRoot = new osg::Group;

   osg::Sphere* sphere = new osg::Sphere(osg::Vec3(t_t[0], t_t[1], t_t[2]), 0.3);
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


      auto R_1_d = exp(-0.5*m_1);
      auto R_2_d = exp(-0.5*m_2);
      auto R_1 = eval(R_1_d*R_1_0);
      auto R_2 = eval(R_2_d*R_2_0);
      double alpha = 0.01;

      //Ball joint R_1:
      auto grad_m1_U = K_1*m_1 + lambda_R + lambda_t*grade<2>(0.5*R_1*t_1*t_1*~R_1 + 0.5*R_1*R_2*t_2*t_2*~R_2*~R_1);

      m_1 = m_1 - alpha*grad_m1_U;

      //Ball joint R_2:
      auto grad_m2_U = K_2*m_2 + lambda_R + lambda_t*grade<2>(0.5*R_1*R_2*t_2*t_2*~R_2*~R_1);
      std::cout << "m_2: " << "lambda: " << lambda_t << ", t_2*t_2: " << t_2*t_2 << ", rest: " << R_1*R_2*t_2*t_2*~R_2*~R_1 << std::endl;

      m_2 = m_2 - alpha*grad_m2_U;

      //Lambda
      auto grad_lambdaR_U = 2.0*grade<2>(log(R_1*R_2*~R_t));
      lambda_R = lambda_R - alpha*grad_lambdaR_U;

      auto grad_lambdat_U = magnitude(grade<1>(R_1*t_1*~R_1 + R_1*R_2*t_2*~R_2*~R_1 - t_t));
      lambda_t = lambda_t - eval(alpha*grad_lambdat_U);

      auto x_1 = eval(grade<1>(R_1*t_1*~R_1));
      auto x_2 = eval(grade<1>(R_1*R_2*t_2*~R_2*~R_1));
      auto R_12 = eval(R_1*R_2);
      //std::cout << "X_1: " << X_1 << ", X_2: " << X_2 << std::endl;
      rod1Transform->setAttitude(osg::Quat(-R_1[3], R_1[2], -R_1[1], R_1[0]));
      rod2Transform->setPosition(osg::Vec3d(x_1[0], x_1[1], x_1[2]));
      rod2Transform->setAttitude(osg::Quat(-R_12[3], R_12[2], -R_12[1], R_12[0]));

      std::cout << "i=" << counter << ": m_1: " << m_1 << ", m_2: " << m_2 << std::endl;
      std::cout << "\tgrad_m1_U: " << grad_m1_U << std::endl;
      std::cout << "\tgrad_m2_U: " << grad_m2_U << std::endl;

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
