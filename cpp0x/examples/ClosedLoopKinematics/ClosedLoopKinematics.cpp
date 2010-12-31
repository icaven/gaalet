///----------------------------------
///Author: Florian Seybold, 2010
///www.hlrs.de
///----------------------------------

#include "cga_osg.h"

#include <osgViewer/Viewer>
#include <osgGA/StateSetManipulator>
#include <osgGA/TrackballManipulator>
#include <osg/PositionAttitudeTransform>
#include <osg/ShapeDrawable>

using namespace cga;

int main()
{
   ::osg::Group* sceneRoot = new ::osg::Group;

   auto x1 = 1.0*e1;
   double r = 0.5;
   cga::osg::Sphere* s1_p = new cga::osg::Sphere(eval(x1 + 0.5*(eval(x1&x1)-r*r)*einf + e0));
   cga::osg::Sphere& s1 = *s1_p;

   ::osg::ShapeDrawable* s1Drawable = new ::osg::ShapeDrawable(s1_p);
   ::osg::Geode* s1Geode = new ::osg::Geode();
   s1Geode->addDrawable(s1Drawable);

   sceneRoot->addChild(s1Geode);


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

   double frameTime = 0.0;
   double sumFrameTime = 0.0;
   double minFrameTime = 0.0;
   double timer = 0.0;
   unsigned int counter = 0;
   while(!viewer.done()) {
      ::osg::Timer_t startFrameTick = ::osg::Timer::instance()->tick();

      s1.assign(x1 + 0.5*(eval(x1&x1)-r*r*timer*timer)*einf + e0);
      s1.update();
      s1Drawable->dirtyDisplayList();

      viewer.frame();

      //work out if we need to force a sleep to hold back the frame rate
      ::osg::Timer_t endFrameTick = ::osg::Timer::instance()->tick();
      frameTime = ::osg::Timer::instance()->delta_s(startFrameTick, endFrameTick);

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
   
   return 0;
}
