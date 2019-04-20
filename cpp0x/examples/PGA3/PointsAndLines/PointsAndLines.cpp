#include "gaalet.h"
#include "pga3.h"

//#include <memory>

#include "OSG_Utilities.h"

// Re-implementation of https://enkimute.github.io/ganja.js/examples/coffeeshop.html#pga3d_points_and_lines
int main()
{

    // Create 5 points and there will be some joining lines between them
    auto origin = pga3::make_point(0, 0, 0);
    auto A = pga3::make_point(0, -1, 0);
    auto B = pga3::make_point(1, 1, -1);
    auto C = pga3::make_point(-1, 1, -1);
    auto D = pga3::make_point(1, 1, 1);
    auto E = pga3::make_point(-1, 1, 1);
//    auto A = pga3::make_point(-1, 0, 0);
//    auto B = pga3::make_point(1, 0, 0);
//    auto C = pga3::make_point(0, -1, 0);
//    auto D = pga3::make_point(0, 1, 0);
//    auto E = pga3::make_point(0, 0, -1);
//    auto F = pga3::make_point(0, 0, 1);

    auto centroid = normalize(A + B + C + D + E);  // Always normalize after operations with the pga3 entities
    
    auto half_way_AB = normalize((A + B) * 0.5);
    auto camera = 0.0 * pga3::e0;

    // Graph the 3D items
    osg::Group* sceneRoot = new osg::Group;
    osg::Geode* cubeGeode = new osg::Geode();
    osg::PositionAttitudeTransform* cubeTransform = new osg::PositionAttitudeTransform();
    cubeTransform->addChild(cubeGeode);
    sceneRoot->addChild(cubeTransform);

    const float point_size = 0.05f;
    std::vector<osg::ShapeDrawable*> drawables;
    drawables.push_back(new_drawable_point(origin, white(), point_size));
    drawables.push_back(new_drawable_point(A, cyan(), point_size));
    drawables.push_back(new_drawable_point(B, red(), point_size));
    drawables.push_back(new_drawable_point(C, magenta(), point_size));
    drawables.push_back(new_drawable_point(D, green(), point_size));
    drawables.push_back(new_drawable_point(E, yellow(), point_size));
    drawables.push_back(new_drawable_point(centroid, grey(0.5), point_size));
//    drawables.push_back(new_drawable_point(F, blue()));
    
    drawables.push_back(new_drawable_arrow(origin, A));
    drawables.push_back(new_drawable_arrow(A, B, cyan(0.25f)));
    drawables.push_back(new_drawable_arrow(A, C));
    drawables.push_back(new_drawable_arrow(A, D));
    drawables.push_back(new_drawable_arrow(B, C));
    drawables.push_back(new_drawable_arrow(B, D, red(0.25f)));
    drawables.push_back(new_drawable_arrow(C, E));
    drawables.push_back(new_drawable_arrow(A, E));
    drawables.push_back(new_drawable_arrow(E, D));

/*
//    std::cout << "Line from origin to A (-x)" << std::endl;
    drawables.push_back(new_drawable_arrow(origin, A, cyan()));
    std::cout << "Line from origin to B (+x)" << std::endl;
    drawables.push_back(new_drawable_arrow(origin, B, red()));
    
//    std::cout << "Line from origin to C (-y)" << std::endl;
    drawables.push_back(new_drawable_arrow(origin, C, magenta()));
    std::cout << "Line from origin to D (+y)" << std::endl;
    drawables.push_back(new_drawable_arrow(origin, D, green()));

    std::cout << "Line from origin to E (-z)" << std::endl;
    drawables.push_back(new_drawable_arrow(origin, E, yellow()));
    std::cout << "Line from origin to F (+z)" << std::endl;
    drawables.push_back(new_drawable_arrow(origin, F, blue()));
*/
//    drawables.push_back(new_drawable_line(F, B, red()));
//    drawables.push_back(new_drawable_line(F, D, green()));

//    auto vertical_line = pga3::line_from_points(origin, F);
//    auto line_O_B = pga3::line_from_points(origin, B);
//    auto translate_O_F = pga3::translator(vertical_line, ::magnitude(pga3::line_from_points(origin, B)));
//    auto moved_O_F = pga3::sandwich(line_O_B, translate_O_B);
    
//    drawables.push_back(new_drawable_line(moved_O_F, , blue()));
    
//    auto translate_B_vertically = pga3::translator(-1 * pga3::k, ::magnitude(pga3::line_from_points(origin, F)));
//    std::cout << "pga3::sandwich(B, translate_B_vertically) " << (pga3::sandwich(B, translate_B_vertically)) << std::endl;
//    drawables.push_back(new_drawable_arrow(B, pga3::sandwich(B, translate_B_vertically), yellow()));
//
//    auto translate_vertically = pga3::translator(-1* pga3::k, ::magnitude(pga3::line_from_points(origin, F)));
//    std::cout << "pga3::sandwich(B, translate_B_vertically) " << (pga3::sandwich(origin, translate_vertically)) << std::endl;
//    drawables.push_back(new_drawable_arrow(pga3::sandwich(origin, translate_vertically), pga3::sandwich(B, translate_vertically), red()));
//
/*
    std::vector<osg::ShapeDrawable*> drawables;
    std::cout << "point from line and plane: " << std::endl;
    auto plane3 = pga3::plane_from_points(A, D, B);

    drawables.push_back(new_drawable_plane(A, D, B));

    
    auto line1 = pga3::line_from_points(centroid, C+E);
    std::cout << "line from points: " << line1 << std::endl;
    drawables.push_back(new_drawable_line(centroid, C+E));

    auto plane1 = pga3::plane_from_points(A,B,C);
    auto plane2 = pga3::plane_from_points(B,C,D);
    std::cout << "line from planes: plane1: " << plane1 << " plane2: " << plane2 << std::endl;
//    drawables.push_back(new_drawable_line(pga3::line_from_planes(plane1, plane2), cyan(0.025f)));
    
    auto P = pga3::point_from_line_and_plane(line1, plane3);
    std::cout << "point_from_line_and_plane " << P << std::endl;
 */   
    for(auto d : drawables) {
        cubeGeode->addDrawable(d);
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
