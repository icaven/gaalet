
#include "OSG_Utilities.h"

inline double sq(double x)
{
    return x * x;
}

// Re-implementation of https://enkimute.github.io/ganja.js/examples/coffeeshop.html#pga3d_icosahedron
int main()
{
    // Our camera position and orientation
    auto camera = 0.0 * pga3::e0;

    // We construct faces, edges and vertices of an icosahedron.
    auto r = pga3::rotor(pga3::e13, M_PI / 2.5);
    pga3::Point_t A = pga3::make_point(0., 1., 0.);
    pga3::Point_t B = pga3::make_point(sqrt(1. - sq(atan(0.5))), atan(0.5), 0.);
    pga3::Point_t C = pga3::sandwich(pga3::sandwich(B, pga3::e2), pga3::rotor(pga3::e13, M_PI / 5));
    pga3::Point_t D = pga3::sandwich(A, pga3::e2);

    // Graph the 3D items
    osg::Group* sceneRoot = new osg::Group;
    osg::Geode* cubeGeode = new osg::Geode();
    osg::PositionAttitudeTransform* cubeTransform = new osg::PositionAttitudeTransform();
    cubeTransform->addChild(cubeGeode);
    sceneRoot->addChild(cubeTransform);

    std::vector<osg::ShapeDrawable*> drawables;
    const float VERTEX_RADIUS = 0.02f;
    drawables.push_back(new_drawable_point(A, white(), VERTEX_RADIUS));

    // vertices
    auto vertex_colour = colour(float(0x44) / 255., float(0x44) / 255., 1.0);
    for(int i = 0; i < 5; i++) {
        B = pga3::sandwich(B, r);
        drawables.push_back(new_drawable_point(B, vertex_colour, VERTEX_RADIUS));
        C = pga3::sandwich(C, r);
        drawables.push_back(new_drawable_point(C, vertex_colour, VERTEX_RADIUS));
        drawables.push_back(new_drawable_point(D, vertex_colour, VERTEX_RADIUS));
    }

    // edges
    const float EDGE_THICKNESS = 0.01f;
    auto line_colour = colour(float(0x44) / 255., float(0x44) / 255., float(0x44) / 255.);
    for(int i = 0; i < 5; i++) {
        drawables.push_back(new_drawable_line(A, B, line_colour, EDGE_THICKNESS));
        drawables.push_back(new_drawable_line(B, C, line_colour, EDGE_THICKNESS));
        
        auto next_B = pga3::sandwich(B, r);
        drawables.push_back(new_drawable_line(B, next_B, line_colour, EDGE_THICKNESS));
        B = next_B;

        drawables.push_back(new_drawable_line(B, C, line_colour, EDGE_THICKNESS));

        auto next_C = pga3::sandwich(C, r);
        drawables.push_back(new_drawable_line(C, next_C, line_colour, EDGE_THICKNESS));
        C = next_C;

        drawables.push_back(new_drawable_line(D, C, line_colour, EDGE_THICKNESS));
    }
    
    // faces
    auto face_colour = colour(float(0xff) / 255., float(0xcc) / 255., float(0xcc) / 255.);
    for(int i = 0; i < 5; i++) {
        auto next_B = pga3::sandwich(B, r);
        drawables.push_back(new_drawable_triangle(A, B, next_B, face_colour, true));
        drawables.push_back(new_drawable_triangle(B, next_B, C, face_colour, true));
        B = next_B;
        auto next_C = pga3::sandwich(C, r);
        
        drawables.push_back(new_drawable_triangle(C, B, next_C, face_colour, true));
        drawables.push_back(new_drawable_triangle(C, D, next_C, face_colour, true));
        C = next_C;
    }

    for(auto p : drawables) {
        cubeGeode->addDrawable(p);
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
    while(!viewer.done()) {
        throttle.begin();

        // Update objects and camera position here

        viewer.frame();
        throttle.end();
    }
    return 0;
}
