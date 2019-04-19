#include "gaalet.h"
#include "pga3.h"

//#include <memory>

#include "OSG_Utilities.h"

// Re-implementation of https://enkimute.github.io/ganja.js/examples/coffeeshop.html#pga3d_points_and_lines
int main()
{

    // Create 5 points and there will be some joining lines between them
    auto origin = pga3::point(0, 0, 0);
    auto A = pga3::point(0, -1, 0);
    auto B = pga3::point(1, 1, -1);
    auto C = pga3::point(-1, 1, -1);
    auto D = pga3::point(1, 1, 1);
    auto E = pga3::point(-1, 1, 1);
//    auto A = pga3::point(-1, 0, 0);
//    auto B = pga3::point(1, 0, 0);
//    auto C = pga3::point(0, -1, 0);
//    auto D = pga3::point(0, 1, 0);
//    auto E = pga3::point(0, 0, -1);
//    auto F = pga3::point(0, 0, 1);

    auto centroid = normalize(A + B + C + D + E);  // Always normalize after operations with the pga3 entities
    
    auto half_way_AB = normalize((A + B) * 0.5);
    auto camera = 0.0 * pga3::e0;

    // Graph the 3D items
    osg::Group* sceneRoot = new osg::Group;
    osg::Geode* cubeGeode = new osg::Geode();
    osg::PositionAttitudeTransform* cubeTransform = new osg::PositionAttitudeTransform();
    cubeTransform->addChild(cubeGeode);
    sceneRoot->addChild(cubeTransform);

    std::vector<osg::ShapeDrawable*> drawable_points;
    drawable_points.push_back(new_drawable_point(origin, white()));
    drawable_points.push_back(new_drawable_point(A, cyan()));
    drawable_points.push_back(new_drawable_point(B, red()));
    drawable_points.push_back(new_drawable_point(C, magenta()));
    drawable_points.push_back(new_drawable_point(D, green()));
    drawable_points.push_back(new_drawable_point(E, yellow()));
//    drawable_points.push_back(new_drawable_point(F, blue()));
    drawable_points.push_back(new_drawable_point(centroid, grey(0.5)));
    
    std::vector<osg::ShapeDrawable*> drawable_lines;
    drawable_lines.push_back(new_drawable_line(origin, A));
    drawable_lines.push_back(new_drawable_line(A, B, cyan(0.25f)));
    drawable_lines.push_back(new_drawable_line(A, C));
    drawable_lines.push_back(new_drawable_line(A, D));
    drawable_lines.push_back(new_drawable_line(B, C));
    drawable_lines.push_back(new_drawable_line(B, D, red(0.25f)));
    drawable_lines.push_back(new_drawable_line(C, E));
    drawable_lines.push_back(new_drawable_line(A, E));
    drawable_lines.push_back(new_drawable_line(E, D));

//    std::cout << "Line from origin to A (-x)" << std::endl;
//    drawable_lines.push_back(new_drawable_line(origin, A, cyan()));
//    std::cout << "Line from origin to B (+x)" << std::endl;
//    drawable_lines.push_back(new_drawable_line(origin, B, red()));
//    
//    std::cout << "Line from origin to C (-y)" << std::endl;
//    drawable_lines.push_back(new_drawable_line(origin, C, magenta()));
//    std::cout << "Line from origin to D (+y)" << std::endl;
//    drawable_lines.push_back(new_drawable_line(origin, D, green()));
//
//    std::cout << "Line from origin to E (-z)" << std::endl;
//    drawable_lines.push_back(new_drawable_line(origin, E, yellow()));
//    std::cout << "Line from origin to F (+z)" << std::endl;
//    drawable_lines.push_back(new_drawable_line(origin, F, blue()));
//

    std::vector<osg::ShapeDrawable*> drawable_planes;
    std::cout << "point from line and plane: " << std::endl;
    auto plane3 = pga3::plane_from_points(A, D, B);

    drawable_planes.push_back(new_drawable_plane(A, D, B));

    
    auto line1 = pga3::line_from_points(centroid, C+E);
    std::cout << "line from points: " << line1 << std::endl;
    drawable_lines.push_back(new_drawable_line(centroid, C+E));

    auto plane1 = pga3::plane_from_points(A,B,C);
    auto plane2 = pga3::plane_from_points(B,C,D);
    std::cout << "line from planes: plane1: " << plane1 << " plane2: " << plane2 << std::endl;
//    drawable_lines.push_back(new_drawable_line(line_from_planes(plane1, plane2), cyan(0.025f)));
    
    auto P = pga3::point_from_line_and_plane(line1, plane3);
    std::cout << "point_from_line_and_plane " << P << std::endl;
    
    for(auto p : drawable_points) {
        cubeGeode->addDrawable(p);
    }
    for(auto l : drawable_lines) {
        cubeGeode->addDrawable(l);
    }
    for(auto pl : drawable_planes) {
        cubeGeode->addDrawable(pl);
    }
    // Graph the 3D items
    //    document.body.appendChild(this.graph(()=>{
    //      var time=performance.now()/4000;
    //      camera.set(Math.cos(time)+Math.sin(time)*1e13);                      // rotate around Y
    //      return [0xddaaff,[A,B,C],                                            // graph on face
    //              0xAA88FF,[A,B],[A,C],[A,D],[B,C],[B,D],[C,E],[A,E],[E,D],    // graph all edges
    //              0x444444,A,"A",B,"B",C,"C",D,"D",E,"E",                      // graph all vertices
    //              0xFF8888,C+E,centroid,"sum of points",
    //              0x8888FF,line_from_points(centroid,C+E),"line from points ..",
    //              0x44AA44,line_from_planes(A&B&C,B&C&D),"line from planes ..",
    //              0x4488FF,point_from_line_and_plane(line_from_points(centroid,C+E),A&D&B),"point from line and plane
    //              ..", 0xFFAA66,(B&D)+(C&E),"sum of lines"
    //             ];
    //    },{animate:true,camera}));
    //});

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
    while (!viewer.done()) {
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
