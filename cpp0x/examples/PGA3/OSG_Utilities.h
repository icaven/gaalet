// Utility functions to help draw scenes using OpenSceneGraph

#pragma once

#include <cmath>
#include <float.h>
#include <vector>

#include "gaalet.h"
#include "pga3.h"

#include <osgViewer/Viewer>
#include <osgGA/StateSetManipulator>
#include <osgGA/TrackballManipulator>
#include <osg/PositionAttitudeTransform>
#include <osg/ShapeDrawable>

inline bool isclose(double a, double b, double rtol=1e-05, double atol=DBL_EPSILON) {
    return fabs(a - b) <= (atol + rtol * fabs(b));
}


// Debugging output functions
inline void print_point_info(const pga3::Point_t& p, const std::string& name) {
    std::cout << name << ": " << p 
//              << " polar of " << name << ": " << pga3::polar(p) 
              << " dual of " << name << ": " << pga3::dual(p) 
              << std::endl;
    std::cout << "normalized " << name << ": " << normalize(p) << std::endl;
//    std::cout << name << "[0]: " << p[0] << " " << name << "[1]: " << p[1] << " " 
//              << name << "[2]: " << p[2] << " " << name << "[3]: " << p[3] << std::endl;
}

inline void print_line_info(const pga3::Line_t& l, const std::string& name) {
    std::cout << name << ": " << l 
//              << " polar of " << name << ": " << pga3::polar(l) 
//              << " dual of " << name << ": " << pga3::dual(l) 
              << std::endl;
    std::cout << "normalized " << name << ": " << normalize(l)
////              << "normalized " << name << " squared: " << normalize(l*l) 
              << std::endl;
    std::cout << name << "[0]: " << l[0] << " " << name << "[1]: " << l[1] << " " << name << "[2]: " << l[2] << std::endl;
    std::cout << name << "[3]: " << l[3] << " " << name << "[4]: " << l[4] << " " << name << "[5]: " << l[5] << std::endl;
}

// Visualization of elements and operations

class FrameThrottle
{
public:
    FrameThrottle() : startFrameTick(osg::Timer::instance()->tick()) {};
    
    float time() 
    {
        return startFrameTick;
    };
    
    void begin() 
    {
        startFrameTick = osg::Timer::instance()->tick();
    };
    
    void end()
    {
        // work out if we need to force a sleep to hold back the frame rate
        osg::Timer_t endFrameTick = osg::Timer::instance()->tick();
        frameTime = osg::Timer::instance()->delta_s(startFrameTick, endFrameTick);

        sumFrameTime += frameTime;
        if(counter == 1000) {
//            std::cout << "Average frame time: " << sumFrameTime / 1000.0 << std::endl;
            sumFrameTime = 0.0;
            counter = 0;
        } else {
            counter++;
        }

        timer += frameTime;

        if(frameTime < minFrameTime) {
            OpenThreads::Thread::microSleep(static_cast<unsigned int>(1000000.0 * (minFrameTime - frameTime)));
        }
    };

protected:
    osg::Timer_t startFrameTick;
    double frameTime = 0.0;
    double sumFrameTime = 0.0;
    double minFrameTime = 1.0/60.0; //1.0/120.0; // seconds
    double timer = 0.0;
    unsigned int counter = 0;
};


// An animation loop that uses a callback to allow updates in the loop
class AnimatedScene
{
public:
    AnimatedScene() : sceneRoot(0) {};

    osg::Geode* setup()
    {
        sceneRoot = new osg::Group;
        rootGeode = new osg::Geode();
        worldTransform = new osg::PositionAttitudeTransform();
        worldTransform->addChild(rootGeode);
        sceneRoot->addChild(worldTransform);
        return rootGeode;
    };

    void set_attitude_position(osg::Quat& attitude, osg::Vec3& position) {
        worldTransform->setAttitude(attitude);
        worldTransform->setPosition(position);
        rootGeode->dirtyBound();
    }
    
    osgViewer::Viewer& initialize_viewer() {
        viewer.setSceneData(sceneRoot);
        osg::StateSet* state_set = viewer.getCamera()->getOrCreateStateSet();

        viewer.addEventHandler(new osgGA::StateSetManipulator(state_set));
        if(!viewer.getCameraManipulator() && viewer.getCamera()->getAllowEventFocus()) {
            viewer.setCameraManipulator(new osgGA::TrackballManipulator());
        }
        viewer.setReleaseContextAtEndOfFrameHint(false);

        if(!viewer.isRealized()) {
            viewer.realize();
        }
        return viewer;
    }
    
    void loop(void (*loop_update)(float frame_time, void* data), void* user_data)
    {
        if (!viewer.isRealized()) {
            initialize_viewer();
        }

        // Animation loop:
        FrameThrottle throttle;
        while(!viewer.done()) {
            throttle.begin();

            // Update objects and camera here
            (*loop_update)(throttle.time(), user_data);
            
            viewer.frame();
            throttle.end();
        }
    };

    osgViewer::Viewer viewer;
    osg::Group* sceneRoot;
    osg::Geode* rootGeode;
    osg::PositionAttitudeTransform* worldTransform;
};

// Colours with selectable opacity
inline osg::Vec4 red(float alpha = 1.f)
{
    return osg::Vec4(1, 0, 0, alpha);
}
inline osg::Vec4 green(float alpha = 1.f)
{
    return osg::Vec4(0, 1, 0, alpha);
}
inline osg::Vec4 blue(float alpha = 1.f)
{
    return osg::Vec4(0, 0, 1, alpha);
}
inline osg::Vec4 yellow(float alpha = 1.f)
{
    return osg::Vec4(1, 1, 0, alpha);
}
inline osg::Vec4 magenta(float alpha = 1.f)
{
    return osg::Vec4(1, 0, 1, alpha);
}
inline osg::Vec4 cyan(float alpha = 1.f)
{
    return osg::Vec4(0, 1, 1, alpha);
}
inline osg::Vec4 grey(float level = 1.f, float alpha = 1.f)
{
    return osg::Vec4(level, level, level, alpha);
}

inline osg::Vec4 white(float alpha = 1.f)
{
    return osg::Vec4(1, 1, 1, alpha);
}

inline osg::Vec4 colour(float r, float g, float b, float alpha = 1.f)
{
    return osg::Vec4(r, g, b, alpha);
}

// Function to convert a pga3 point to OSG Vec3
inline osg::Vec3 Vec3(const pga3::Point_t& p)
{
    return osg::Vec3(pga3::Point_x(p), pga3::Point_y(p), pga3::Point_z(p));
}

// Function to convert n angle and pga3 line to OSG Quat for the attitude
inline osg::Quat Quat(double angle, const pga3::Line_t& line)
{
    pga3::Line_t n_line = sin(angle) * normalize(line);
    return osg::Quat(n_line.template element<pga3::I_CONF>(), 
                     n_line.template element<pga3::J_CONF>(),
                     n_line.template element<pga3::K_CONF>(), 
                     cos(angle));
}

// Function to convert a OSG Quat to a PGA3 rotor
inline auto Quat2Rotor(const osg::Quat q)
{
    double angle;
    osg::Vec3d axis;
    q.getRotate(angle, axis);
    auto pga_axis = axis[0] * pga3::i + axis[1] * pga3::j + axis[2] * pga3::k;
    auto R = pga3::rotor(pga_axis, angle);
//    std::cout << "Rotor: " << R << "R*~R: " << R*(~R) << std::endl;
    return R;
}

const float DEFAULT_RADIUS_OF_DRAWN_POINT = 0.1f;
const float DEFAULT_LINE_THICKNESS = 0.0125f;
const float DEFAULT_THICKNESS_OF_PLANE = 0.025f;


// Functions to create drawable primitives

osg::ShapeDrawable* new_drawable_point(const pga3::Point_t& p, const osg::Vec4& colour = white(),
                                       float point_radius=DEFAULT_RADIUS_OF_DRAWN_POINT);

osg::ShapeDrawable* new_drawable_line(const pga3::Point_t& start_pt, const pga3::Point_t& end_pt, 
                                      const osg::Vec4& colour = grey(0.5),
                                      const float line_thickness=DEFAULT_LINE_THICKNESS);
osg::CompositeShape* new_arrow(const pga3::Point_t& start_pt, const pga3::Point_t& end_pt, 
                               const float line_thickness=DEFAULT_LINE_THICKNESS);
osg::ShapeDrawable* new_drawable_arrow(const pga3::Point_t& start_pt, const pga3::Point_t& end_pt, 
                                       const osg::Vec4& colour = grey(0.5),
                                       const float line_thickness=DEFAULT_LINE_THICKNESS);
osg::ShapeDrawable* new_drawable_arrow(const pga3::Point_t& start_pt, const pga3::Line_t& line, 
                                       const osg::Vec4& colour = grey(0.5),
                                       const float line_thickness=DEFAULT_LINE_THICKNESS);

osg::ShapeDrawable* new_drawable_plane(const pga3::Point_t& p1, const pga3::Point_t& p2, const pga3::Point_t& p3,
                                       const osg::Vec4& colour = grey(0.5), const float plane_thickness=DEFAULT_THICKNESS_OF_PLANE);

osg::ShapeDrawable* new_drawable_triangle(const pga3::Point_t& p1, const pga3::Point_t& p2, const pga3::Point_t& p3,
                                          const osg::Vec4& colour = grey(0.5), bool draw_normal=false);
                                       
                                       
