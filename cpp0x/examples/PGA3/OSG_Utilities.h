// Utility functions to help draw scenes using OpenSceneGraph

#include <cmath>
#include <vector>

#include "gaalet.h"
#include "pga3.h"

#include <osgViewer/Viewer>
#include <osgGA/StateSetManipulator>
#include <osgGA/TrackballManipulator>
#include <osg/PositionAttitudeTransform>
#include <osg/ShapeDrawable>

bool isclose(double a, double b, double rtol=1e-05, double atol=DBL_EPSILON) {
    return abs(a - b) <= (atol + rtol * abs(b));
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
            std::cout << "Average frame time: " << sumFrameTime / 1000.0 << std::endl;
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
    std::cout << "Rotor: " << R << "R*~R: " << R*(~R) << std::endl;
    return R;
}


const float DEFAULT_RADIUS_OF_DRAWN_POINT = 0.1f;

// Functions to create drawable primitives
osg::ShapeDrawable* new_drawable_point(const pga3::Point_t& p, const osg::Vec4& colour = white(),
                                       float point_radius=DEFAULT_RADIUS_OF_DRAWN_POINT)
{
    osg::Sphere* sphere = new osg::Sphere(Vec3(p), point_radius);
    osg::ShapeDrawable* drawable = new osg::ShapeDrawable(sphere);
    drawable->setColor(colour);
    return drawable;
}


// Debugging output functions
void print_point_info(const pga3::Point_t& p, const std::string& name) {
    std::cout << name << ": " << p 
//              << " polar of " << name << ": " << pga3::polar(p) 
              << " dual of " << name << ": " << pga3::dual(p) 
              << std::endl;
    std::cout << "normalized " << name << ": " << normalize(p) << std::endl;
//    std::cout << name << "[0]: " << p[0] << " " << name << "[1]: " << p[1] << " " 
//              << name << "[2]: " << p[2] << " " << name << "[3]: " << p[3] << std::endl;
}

void print_line_info(const pga3::Line_t& l, const std::string& name) {
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

const float DEFAULT_LINE_THICKNESS = 0.0125f;


osg::ShapeDrawable*
new_drawable_line(const pga3::Point_t& start_pt, const pga3::Point_t& end_pt, 
                  const osg::Vec4& colour = grey(0.5),
                  const float line_thickness=DEFAULT_LINE_THICKNESS)
{
    pga3::Line_t line = pga3::line_from_points(start_pt, end_pt);
    osg::Vec3 direction = Vec3(end_pt) - Vec3(start_pt);
    float length = direction.length();
              
    // OSG Cylinders always start out in the vertical direction (+Z-axis) and centered on the starting point
    // so rotate and translate as needed
    osg::Cylinder* cylinder = new osg::Cylinder(Vec3(start_pt), line_thickness, length);

    // Determine the angle between the line described by the start and end points and the Z-axis
    double angle = acos((pga3::k & normalize(line)).template element<0>());
    if (isclose(angle, 0.)) {
        // Parallel to the Z-axis, just shift along that axis in the positive direction
        cylinder->setCenter(cylinder->getCenter() + osg::Vec3(0., 0., length*0.5f));
    }
    else if (isclose(angle, M_PI)) {
        // Parallel to the Z-axis, just shift along that axis in the negative direction
        cylinder->setCenter(cylinder->getCenter() + osg::Vec3(0., 0., -length*0.5f));
    }
    else {
        // Determine the plane that the Z-axis and the line forms and then compute the perpendicular
        // to that plane; the angle will correspond to the angle around that perpendicular
        auto end_of_cylinder = normalize(pga3::sandwich(start_pt, pga3::translator(pga3::k, length)));
        auto plane_with_lines = pga3::plane_from_points(start_pt,  end_pt, end_of_cylinder);
        auto perpendicular_to_plane = start_pt & plane_with_lines;
        
        // Rotate the cylinder to be in the direction of the line
        osg::Quat q = Quat(angle/2.0, perpendicular_to_plane);
        cylinder->setRotation(q);
    
        // Translate the center point
        auto t = pga3::translator(line, -length*0.5);
        auto new_center = pga3::sandwich(start_pt, t);
        cylinder->setCenter(Vec3(new_center));
    }

    osg::ShapeDrawable* drawable = new osg::ShapeDrawable(cylinder);
    drawable->setColor(colour);
    return drawable;
}

osg::ShapeDrawable*
new_drawable_arrow(const pga3::Point_t& start_pt, const pga3::Point_t& end_pt, 
                  const osg::Vec4& colour = grey(0.5),
                  const float line_thickness=DEFAULT_LINE_THICKNESS)
{
    pga3::Line_t line = pga3::line_from_points(start_pt, end_pt);
    osg::Vec3 origin = Vec3(start_pt);
    osg::Vec3 direction = Vec3(end_pt) - origin;
    float length = direction.length();
              
    // OSG Cylinders and Cones always start out in the vertical direction (+Z-axis) 
    // and origined on the starting point so rotate and translate as needed
    float arrow_head_length = 5.f*line_thickness;
    float arrow_base_radius = 3.f*line_thickness;
    osg::Cone* arrow_head = new osg::Cone(origin, arrow_base_radius, 5.f*line_thickness);

    // Determine the offset of the tip from the location of the cone origin
    float head_offset = (1.f - arrow_head->getBaseOffsetFactor()) * arrow_head_length;
    float shaft_length = length - head_offset;
    osg::Cylinder* shaft = new osg::Cylinder(origin, line_thickness, shaft_length);

    // Determine the angle between the line described by the start and end points and the Z-axis
    double angle = acos((pga3::k & normalize(line)).template element<0>());
    if (isclose(angle, 0.)) {
        // Parallel to the Z-axis, just shift along that axis in the positive direction
        shaft->setCenter(shaft->getCenter() + osg::Vec3(0., 0., shaft_length*0.5f));
        arrow_head->setCenter(shaft->getCenter() + 
            osg::Vec3(0., 0., fmax((length-head_offset)*0.5f, arrow_head_length)));
    }
    else if (isclose(angle, M_PI)) {
        // Parallel to the Z-axis, just shift along that axis in the negative direction
        shaft->setCenter(shaft->getCenter() + osg::Vec3(0., 0., -shaft_length*0.5f));
        arrow_head->setRotation(osg::Quat(M_PI, osg::Vec3d(1., 0., 0.))); // Reflect the arrow head
        arrow_head->setCenter(shaft->getCenter() + 
            osg::Vec3(0., 0., fmin(-(length-head_offset)*0.5f, -arrow_head_length)));
   }
    else {
        // Determine the plane that the Z-axis and the line forms and then compute the perpendicular
        // to that plane; the angle will correspond to the angle around that perpendicular
        auto end_of_shaft = normalize(pga3::sandwich(start_pt, pga3::translator(pga3::k, shaft_length)));
        auto plane_with_lines = pga3::plane_from_points(start_pt,  end_pt, end_of_shaft);
        auto perpendicular_to_plane = start_pt & plane_with_lines;
        
        // Rotate the shaft to be in the direction of the line
        osg::Quat q = Quat(angle/2.0, perpendicular_to_plane);
        shaft->setRotation(q);

        // Translate the center point of the arrow shaft
        auto t = pga3::translator(line, -shaft_length*0.5);
        auto new_center = pga3::sandwich(start_pt, t);
        shaft->setCenter(Vec3(new_center));
        
        // Rotate and translate the arrow head
        arrow_head->setRotation(q);
        auto arrow_translator = pga3::translator(line, (double) -fmax(length-head_offset, arrow_head_length));
        auto arrow_center = pga3::sandwich(start_pt, arrow_translator);
        arrow_head->setCenter(Vec3(arrow_center));
    }

    osg::CompositeShape* arrow = new osg::CompositeShape();
    arrow->addChild(shaft);
    arrow->addChild(arrow_head);
    
    osg::ShapeDrawable* drawable = new osg::ShapeDrawable(arrow);
    drawable->setColor(colour);
    return drawable;
}

const float DEFAULT_THICKNESS_OF_PLANE = 0.025f;

osg::ShapeDrawable*
new_drawable_plane(const pga3::Point_t& p1, const pga3::Point_t& p2, const pga3::Point_t& p3,
                   const osg::Vec4& colour = grey(0.5), const float plane_thickness=DEFAULT_THICKNESS_OF_PLANE)
{
    pga3::Plane_t the_plane = pga3::plane_from_points(p1, p2, p3);
    std::cout << "The plane from points: " << the_plane << std::endl;
    osg::Plane implicit_plane = osg::Plane(the_plane[1], the_plane[2], the_plane[3], the_plane[0]);
    osg::Vec4 v = implicit_plane.asVec4();
    std::cout << "The plane implicitly from the points: " << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << std::endl;

    osg::Plane plane_directly = osg::Plane(Vec3(p1), Vec3(p2), Vec3(p3));;
    v = plane_directly.asVec4();
    std::cout << "The plane directly from the points: " << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << std::endl;

    double x_length = eval(::magnitude(pga3::line_from_points(p1, p2)));
    double y_length = eval(::magnitude(pga3::line_from_points(p1, p3)));
    double angle = acos((pga3::i & normalize(pga3::line_from_points(p1, p2))).template element<0>());
    std::cout << angle << " " << x_length << " " << y_length << std::endl;

    osg::Box* box_as_plane = new osg::Box(Vec3(p1),float(x_length),float(y_length), plane_thickness);
    auto perpendicular_to_plane = p1 & the_plane;

    osg::Quat q = Quat(angle/2.0, perpendicular_to_plane);
    box_as_plane->setRotation(q);
    osg::ShapeDrawable* drawable = new osg::ShapeDrawable(box_as_plane);
    
    
//    drawable->setColor(colour);
    drawable->setColor(osg::Vec4(0.5, 0.5, 0.5, 0.1));
    return drawable;
}