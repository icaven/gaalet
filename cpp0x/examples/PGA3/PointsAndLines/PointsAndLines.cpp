
#include "OSG_Utilities.h"

// Re-implementation of https://enkimute.github.io/ganja.js/examples/coffeeshop.html#pga3d_points_and_lines
int main()
{
    // Create 5 points and there will be some joining lines between them
    pga3::Point_t origin = pga3::make_point(0, 0, 0);
    pga3::Point_t point_on_x = pga3::make_point(1, 0, 0);
    pga3::Point_t point_on_y = pga3::make_point(0, 1, 0);
    pga3::Point_t point_on_z = pga3::make_point(0, 0, 1);
    pga3::Point_t A = pga3::make_point(0, -1, 0);
    pga3::Point_t B = pga3::make_point(1, 1, -1);
    pga3::Point_t C = pga3::make_point(-1, 1, -1);
    pga3::Point_t D = pga3::make_point(1, 1, 1);
    pga3::Point_t E = pga3::make_point(-1, 1, 1);

    // Need to normalize when computing the centroid
    auto centroid = pga3::normalize(A + B + C + D + E);
    
//    auto camera = 0.0 * pga3::e0;

    // Graph the 3D items
    auto* sceneRoot = new osg::Group;
    auto* cubeGeode = new osg::Geode();
    auto* cubeTransform = new osg::PositionAttitudeTransform();
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

    // Draw in the x, y, and z axes
    drawables.push_back(new_drawable_arrow(origin, pga3::Point_t(origin + pga3::E1), red()));
    drawables.push_back(new_drawable_arrow(origin, pga3::Point_t(origin + pga3::E2), green()));
//    drawables.push_back(new_drawable_arrow(origin, pga3::Point_t(origin + pga3::E3), blue()));

    drawables.push_back(new_drawable_point(pga3::Point_t(origin + pga3::E1), red(), point_size));
    drawables.push_back(new_drawable_point(pga3::Point_t(origin + pga3::E2), green(), point_size));
    drawables.push_back(new_drawable_point(pga3::Point_t(origin + pga3::E3), blue(), point_size));

    drawables.push_back(new_drawable_arrow(origin, pga3::Point_t(origin + pga3::e1 * pga3::I), red()));
    drawables.push_back(new_drawable_arrow(origin, pga3::Point_t(origin + pga3::e2 * pga3::I), green()));
    drawables.push_back(new_drawable_arrow(origin, pga3::Point_t(origin + pga3::e3 * pga3::I), blue()));

    drawables.push_back(new_drawable_triangle(origin, origin+pga3::E2, origin+pga3::E1, red(0.5), true));
    drawables.push_back(new_drawable_triangle(origin, origin+pga3::E1, origin+pga3::E3, green(0.5), true));
    drawables.push_back(new_drawable_triangle(origin, origin+pga3::E3, origin+pga3::E2, blue(0.5), true));

    drawables.push_back(new_drawable_arrow(origin, A));
    drawables.push_back(new_drawable_arrow(A, B, cyan(0.25f)));
    drawables.push_back(new_drawable_arrow(A, C, green()));
    drawables.push_back(new_drawable_arrow(A, D));
    drawables.push_back(new_drawable_arrow(B, C));
    drawables.push_back(new_drawable_arrow(B, D, red(0.25f)));
    drawables.push_back(new_drawable_arrow(C, E));
    drawables.push_back(new_drawable_arrow(A, E));
    drawables.push_back(new_drawable_arrow(E, D));

//    drawables.push_back(new_drawable_plane(A, B, D));
    drawables.push_back(new_drawable_triangle(A, B, D, grey(0.5), true));

    for(auto d : drawables) {
        cubeGeode->addDrawable(d);
    }

    // Ensure that the window is displayed on only one monitor (instead of being split across them)
    osgViewer::Viewer viewer;
    osgViewer::ViewerBase::Views views;
    viewer.getViews(views);
    osg::ref_ptr<osgViewer::SingleWindow> win = new osgViewer::SingleWindow(20,30, 1800, 1000, 0);
    views[0]->apply(win);

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
