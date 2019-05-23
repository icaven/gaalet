#include "gaalet.h"
#include "pga3.h"

//#include <memory>

#include "OSG_Utilities.h"

// Re-implementation of https://enkimute.github.io/ganja.js/examples/coffeeshop.html#pga3d_points_and_lines
int main()
{
    // Create 5 points and there will be some joining lines between them
    pga3::Point_t origin = pga3::make_point(0, 0, 0);
    pga3::Point_t A = pga3::make_point(0, -1, 0);
    pga3::Point_t B = pga3::make_point(1, 1, -1);
    pga3::Point_t C = pga3::make_point(-1, 1, -1);
    pga3::Point_t D = pga3::make_point(1, 1, 1);
    pga3::Point_t E = pga3::make_point(-1, 1, 1);
    
    // Need to normalize when computing the centroid
    auto centroid = normalize(A + B + C + D + E);  
    
//    auto camera = 0.0 * pga3::e0;

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
    
    // Use arrows to show the direction of the lines
    drawables.push_back(new_drawable_arrow(origin, A));
    drawables.push_back(new_drawable_arrow(A, B, cyan(0.25f)));
    drawables.push_back(new_drawable_arrow(A, C));
    drawables.push_back(new_drawable_arrow(A, D));
    drawables.push_back(new_drawable_arrow(B, C));
    drawables.push_back(new_drawable_arrow(B, D, red(0.25f)));
    drawables.push_back(new_drawable_arrow(C, E));
    drawables.push_back(new_drawable_arrow(A, E));
    drawables.push_back(new_drawable_arrow(E, D));

//    drawables.push_back(new_drawable_plane(A, B, D));
    drawables.push_back(new_drawable_triangle(A, B, D, grey(0.5), true));
    
    auto sum_of_lines = pga3::line_from_points(B, C) + pga3::line_from_points(C, E);
    std::cout << "line_from_points(B, C): " << pga3::line_from_points(B, C) << std::endl;
    std::cout << "line_from_points(C, E): " << pga3::line_from_points(C, E) << std::endl;
    std::cout << "sum of lines: " << sum_of_lines << std::endl;
    std::cout << "sum of lines normalized: " << ::normalize(sum_of_lines) << std::endl;

    for(auto d : drawables) {
        cubeGeode->addDrawable(d);
    }

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
        
        // Update objects and camera here
        
        viewer.frame();
        throttle.end();
    }
    return 0;
}
