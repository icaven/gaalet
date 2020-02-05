// Utility functions to help draw scenes using OpenSceneGraph

#pragma once

#include <cmath>
#include <cfloat>
#include <vector>

#include "gaalet.h"
#include "pga3.h"

#include <osgViewer/Viewer>
#include <osgViewer/config/SingleWindow>
#include <osgViewer/CompositeViewer>
#include <osgGA/StateSetManipulator>
#include <osgGA/TrackballManipulator>
#include <osg/PositionAttitudeTransform>
#include <osg/ShapeDrawable>

inline bool isclose(double a, double b, double rtol=1e-05, double atol=DBL_EPSILON) {
    return fabs(a - b) <= (atol + rtol * fabs(b));
}

// Visualization of elements and operations

class FrameThrottle
{
public:
    FrameThrottle() : m_startFrameTick(osg::Timer::instance()->tick()) {};
    
    float time() 
    {
        return m_startFrameTick;
    };
    
    void begin() 
    {
        m_startFrameTick = osg::Timer::instance()->tick();
    };
    
    void end()
    {
        // work out if we need to force a sleep to hold back the frame rate
        osg::Timer_t endFrameTick = osg::Timer::instance()->tick();
        m_frameTime = osg::Timer::instance()->delta_s(m_startFrameTick, endFrameTick);

        m_sumFrameTime += m_frameTime;
        if(m_counter == 1000) {
//            std::cout << "Average frame time: " << sumFrameTime / 1000.0 << std::endl;
            m_sumFrameTime = 0.0;
            m_counter = 0;
        } else {
            m_counter++;
        }

        m_timer += m_frameTime;

        if(m_frameTime < m_minFrameTime) {
            OpenThreads::Thread::microSleep(static_cast<unsigned int>(1000000.0 * (m_minFrameTime - m_frameTime)));
        }
    };

protected:
    osg::Timer_t m_startFrameTick;
    double m_frameTime = 0.0;
    double m_sumFrameTime = 0.0;
    double m_minFrameTime = 1.0 / 60.0; //1.0/120.0; // seconds
    double m_timer = 0.0;
    unsigned int m_counter = 0;
};


// An animation loop that uses a callback to allow updates in the loop
class AnimatedScene
{
public:
    AnimatedScene() : m_sceneRoot(nullptr) {};

    osg::Geode* setup()
    {
        m_sceneRoot = new osg::Group();
        m_rootGeode = new osg::Geode();
        m_worldTransform = new osg::PositionAttitudeTransform();
        m_worldTransform->addChild(m_rootGeode);
        m_sceneRoot->addChild(m_worldTransform);
        return m_rootGeode;
    };

    void set_attitude_position(osg::Quat& attitude, osg::Vec3& position) {
        m_worldTransform->setAttitude(attitude);
        m_worldTransform->setPosition(position);
        m_rootGeode->dirtyBound();
    }
    
    void initialize_viewer(int x=20, int y=30, int width=1800, int height=1000, unsigned int screen_number=0) {

        // Ensure that the window is displayed on only one monitor (instead of being split across them)
        osgViewer::ViewerBase::Views views;
        m_viewer.getViews(views);
        osg::ref_ptr<osgViewer::SingleWindow> win = new osgViewer::SingleWindow(x, y, width, height, screen_number);
        views[0]->apply(win);

        m_viewer.setSceneData(m_sceneRoot);
        osg::StateSet* state_set = m_viewer.getCamera()->getOrCreateStateSet();

        m_viewer.addEventHandler(new osgGA::StateSetManipulator(state_set));
        if(!m_viewer.getCameraManipulator() && m_viewer.getCamera()->getAllowEventFocus()) {
            m_viewer.setCameraManipulator(new osgGA::TrackballManipulator());
        }
        m_viewer.setReleaseContextAtEndOfFrameHint(false);

        if(!m_viewer.isRealized()) {
            m_viewer.realize();
        }
    }
    
    void loop(void (*loop_update)(float frame_time, void* data), void* user_data)
    {
        if (!m_viewer.isRealized()) {
            initialize_viewer();
        }

        // Animation loop:
        FrameThrottle throttle;
        while(!m_viewer.done()) {
            throttle.begin();

            // Update objects and camera here
            (*loop_update)(throttle.time(), user_data);

            m_viewer.frame();
            throttle.end();
        }
    };

    osgViewer::Viewer m_viewer;
    osg::Group* m_sceneRoot = nullptr;
    osg::Geode* m_rootGeode = nullptr;
    osg::PositionAttitudeTransform* m_worldTransform = nullptr;
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
    return osg::Quat(n_line.template element<pga3::K_CONF>(),  // x-axis
                     -n_line.template element<pga3::J_CONF>(), // y-axis; need to negate because value is stored negated
                     n_line.template element<pga3::I_CONF>(),  // z-axis
                     cos(angle));
}

// Function to convert a OSG Quat to a PGA3 rotor
inline auto Quat2Rotor(const osg::Quat q)
{
    double angle;
    osg::Vec3d axis;
    q.getRotate(angle, axis);
    auto pga_axis = axis[0] * pga3::k + axis[1] * pga3::j + axis[2] * pga3::i;
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
                                      float line_thickness=DEFAULT_LINE_THICKNESS);
osg::CompositeShape* new_arrow(const pga3::Point_t& start_pt, const pga3::Point_t& end_pt, 
                               float line_thickness=DEFAULT_LINE_THICKNESS);
osg::ShapeDrawable* new_drawable_arrow(const pga3::Point_t& start_pt, const pga3::Point_t& end_pt, 
                                       const osg::Vec4& colour = grey(0.5),
                                       float line_thickness=DEFAULT_LINE_THICKNESS);
osg::ShapeDrawable* new_drawable_arrow(const pga3::Point_t& start_pt, const pga3::Line_t& line, 
                                       const osg::Vec4& colour = grey(0.5),
                                       float line_thickness=DEFAULT_LINE_THICKNESS);

osg::ShapeDrawable* new_drawable_plane(const pga3::Point_t& p1, const pga3::Point_t& p2, const pga3::Point_t& p3,
                                       const osg::Vec4& colour = grey(0.5), float plane_thickness=DEFAULT_THICKNESS_OF_PLANE);

osg::ShapeDrawable* new_drawable_triangle(const pga3::Point_t& p1, const pga3::Point_t& p2, const pga3::Point_t& p3,
                                          const osg::Vec4& colour = grey(0.5), bool draw_normal=false);
                                       
                                       
