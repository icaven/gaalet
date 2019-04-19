
#include "OSG_Utilities.h"

inline double sq(double x) { return x * x; }

// Re-implementation of https://enkimute.github.io/ganja.js/examples/coffeeshop.html#pga3d_icosahedron
int main()
{
    // Our camera position and orientation
    auto camera = 0.0 * pga3::e0;

    // We construct faces, edges and vertices of an icosahedron.
    auto r = pga3::rotor(pga3::e13, M_PI / 2.5);
    pga3::Point_t A = pga3::point(0., 1., 0.);
    pga3::Point_t B = pga3::point(sqrt(1. - sq(atan(0.5))), atan(0.5), 0.);
    pga3::Point_t C = pga3::sandwich(pga3::rotor(pga3::e13, M_PI / 5), pga3::sandwich(pga3::e2, B));
    pga3::Point_t D = pga3::sandwich(pga3::e2, A);

    // Graph the 3D items
    osg::Group* sceneRoot = new osg::Group;
    osg::Geode* cubeGeode = new osg::Geode();
    osg::PositionAttitudeTransform* cubeTransform = new osg::PositionAttitudeTransform();
    cubeTransform->addChild(cubeGeode);
    sceneRoot->addChild(cubeTransform);

    std::vector<osg::ShapeDrawable*> drawable_points;
    drawable_points.push_back(new_drawable_point(A, white()));
    //  drawable_points.push_back(new_drawable_point(B, white()));
    //  drawable_points.push_back(new_drawable_point(C, white()));

    // vertices
    // items.push(0x4444FF);
    auto vertex_colour = colour(float(0x44) / 255., float(0x44) / 255., 1.0);
    for(int i = 0; i < 5; i++) {
        //    drawable_points.push_back(new_drawable_point(A, vertex_colour));
        pga3::Point_t newB(pga3::sandwich(r, B));
        B = newB;
        drawable_points.push_back(new_drawable_point(newB, vertex_colour));
        C = pga3::sandwich(r, C);
        drawable_points.push_back(new_drawable_point(C, vertex_colour));
        drawable_points.push_back(new_drawable_point(D, vertex_colour));
    }

    for(auto p : drawable_points) {
        cubeGeode->addDrawable(p);
    }

    return 0; // XXXX DEBUGGING ONLY
              //
              //  // edges
              //  i//  items.push(0x444444);
              //  for (var i=0;i<5;i++) items.push([A,B],[B,C],[B,B=r>>>  return 0; // XXXX DEBUGGING ONLY
  //
  //  // edges
  //  i//  items.push(0x444444);
//  for (var i=0;i<5;i++) items.push([A,B],[B,C],[B,B=r>>>B],[B,C],[C,C=r>>>C],[1e2>>>A,C]);
//  
//  // faces
//  items.push(0xFFCCCC);
//  for (var i=0;i<5;i++) items.push([A,B,r>>>B],[B,B=r>>>B,C],[C,B,r>>>C],[C,1e2>>>A,C=r>>>C]);
//  
//  // Graph the 3D items
//  document.body.appendChild(this.graph(()=>{
//    var time=performance.now()/4000;    
//    camera.set(rotor(1e13,time)*rotor(1e12,time*1.23131));                // animate camera
//    return items.slice(0,1+((Math.floor(time*50))%(items.length+20)));    // show more and more elements
//  },{gl:true,animate:true,camera})); 

    osgViewer::Viewer viewer;
    viewer.setSceneData(sceneRoot);
    viewer.addEventHandler(new osgGA::StateSetManipulator(viewer.getCamera()->getOrCreateStateSet()));
    if(!viewer.getCameraManipulator() && viewer.getCamera()->getAllowEventFocus()) {
        viewer.setCameraManipulator(new osgGA::TrackballManipulator());
    }
    viewer.setReleaseContextAtEndOfFrameHint(false);

    if(!viewer.isRealized()) {
        viewer.realize();
    }

    // Animation loop:
    FrameThrottle throttle;
    while (false && !viewer.done()) {
        throttle.begin();
        
        // Updating new position of cube
        //      auto p_m = eval(grade<1>(D*e0*(~D)));
        //      cubeTransform->setPosition(osg::Vec3(p_m[0], p_m[1], p_m[2]));
        //      cubeTransform->setAttitude(osg::Quat(-D[3], D[2], -D[1], D[0]));

        viewer.frame();
        throttle.end();
    }
    return 0;  
}
