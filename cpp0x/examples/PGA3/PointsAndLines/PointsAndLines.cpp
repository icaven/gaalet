#include "gaalet.h"
#include "pga3.h"

#include <memory>

#include <osgViewer/Viewer>
#include <osgGA/StateSetManipulator>
#include <osgGA/TrackballManipulator>
#include <osg/PositionAttitudeTransform>
#include <osg/ShapeDrawable>

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
    double minFrameTime = 1.0/120.0; // seconds
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

// Function to convert a pga3 point to OSG Vec3
inline osg::Vec3 Vec3(const pga3::Point_t& p)
{
    return osg::Vec3(pga3::Point_x(p), pga3::Point_y(p), pga3::Point_z(p));
}

// Function to convert a pga3 line to OSG Quat for the attitude
inline osg::Quat Quat(const pga3::Line_t& line)
{
    pga3::Line_t n_line = normalize(line);
//    pga3::rotor()
    return osg::Quat(n_line.template element<pga3::J_CONF>(), 
                     n_line.template element<pga3::I_CONF>(), 
                     n_line.template element<pga3::K_CONF>(), 
                     1.0);
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


// Functions to create drawable primitives
osg::ShapeDrawable* new_drawable_point(const pga3::Point_t& p, const osg::Vec4& colour = white())
{
    osg::Sphere* sphere = new osg::Sphere(Vec3(p), 0.1f);
    osg::ShapeDrawable* drawable = new osg::ShapeDrawable(sphere);
    drawable->setColor(colour);
    return drawable;
}

void print_point_info(const pga3::Point_t& p, const std::string& name) {
//    std::cout << name << ": " << p 
//              << " polar of " << name << ": " << pga3::polar(p) 
//              << " dual of " << name << ": " << ::dual(p) 
//              << std::endl;
//    std::cout << "normalized " << name << ": " << normalize(p) << std::endl;
//    std::cout << name << "[0]: " << p[0] << " " << name << "[1]: " << p[1] << " " 
//              << name << "[2]: " << p[2] << " " << name << "[3]: " << p[3] << std::endl;
}

void print_line_info(const pga3::Line_t& l, const std::string& name) {
//    std::cout << name << ": " << l 
//              << " polar of " << name << ": " << pga3::polar(l) 
//              << " dual of " << name << ": " << ::dual(l) 
//              << std::endl;
//    std::cout << "normalized " << name << ": " << normalize(l)
////              << "normalized " << name << " squared: " << normalize(l*l) 
//              << std::endl;
//    std::cout << name << "[0]: " << l[0] << " " << name << "[1]: " << l[1] << " " << name << "[2]: " << l[2] << std::endl;
//    std::cout << name << "[3]: " << l[3] << " " << name << "[4]: " << l[4] << " " << name << "[5]: " << l[5] << std::endl;
}

osg::ShapeDrawable*
new_drawable_line(const pga3::Point_t& start_pt, const pga3::Point_t& end_pt, 
                  const osg::Vec4& colour = grey(0.5),
                  const float line_thickness=0.025f)
{
    print_point_info(start_pt, "start");
    print_point_info(end_pt, "end");
    pga3::Line_t line = pga3::line_from_points(start_pt, end_pt);
    print_line_info(line, "line");
    osg::Vec3 origin = Vec3(start_pt);
    osg::Vec3 direction = Vec3(end_pt) - origin;
    float length = direction.length();
    std::cout << "direction: (" << direction.x() << ", " 
              << direction.y() << ", " << direction.z() << ") length: " << length << std::endl;
    osg::Cylinder* cylinder = new osg::Cylinder(origin, line_thickness, length);
    
//    osg::Quat q = Quat(line);
//    double angle_radians;
//    double x, y, z;
//    q.getRotate(angle_radians, x, y, z);
//    std::cout << angle_radians << " " << x << " " << y << " " << z << std::endl;
//    cylinder->setRotation(q);
//    cylinder->setCenter(cylinder->getCenter() + 
//                        Vec3(pga3::sandwich(pga3::point(0., 0., length * 0.5), Quat2Rotor(q))));
    osg::ShapeDrawable* drawable = new osg::ShapeDrawable(cylinder);
    drawable->setColor(colour);
    return drawable;
}

osg::ShapeDrawable*
new_drawable_plane(const pga3::Point_t& p1, const pga3::Point_t& p2, const pga3::Point_t& p3,
                   const osg::Vec4& colour = grey(0.5), const float plane_thickness=0.025f)
{
    pga3::Plane_t the_plane = pga3::plane_from_points(p1, p2, p3);
    std::cout << "The plane from points: " << the_plane << std::endl;
    osg::Plane implicit_plane = osg::Plane(the_plane[1], the_plane[2], the_plane[3], the_plane[0]);
    osg::Vec4 v = implicit_plane.asVec4();
    std::cout << "The plane implicitly from the points: " << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << std::endl;

    osg::Plane plane_directly = osg::Plane(Vec3(p1), Vec3(p2), Vec3(p3));;
    v = plane_directly.asVec4();
    std::cout << "The plane directly from the points: " << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << std::endl;
    osg::ShapeDrawable* drawable = new osg::ShapeDrawable(new osg::Box());
//    drawable->setColor(colour);
    drawable->setColor(osg::Vec4(0.5, 0.5, 0.5, 0.1));
    return drawable;
}

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
    for(auto p : drawable_points) {
        cubeGeode->addDrawable(p);
    }
    
    std::vector<osg::ShapeDrawable*> drawable_lines;
//    std::cout << "Line from origin to A" << std::endl;
//    drawable_lines.push_back(new_drawable_line(origin, A));
    std::cout << "Line from A to B, cyan" << std::endl;
    drawable_lines.push_back(new_drawable_line(A, B, cyan(0.25f)));
////    drawable_lines.push_back(new_drawable_line(A, C));
////    drawable_lines.push_back(new_drawable_line(A, D));
////    drawable_lines.push_back(new_drawable_line(B, C));
    std::cout << "Line from B to D, red" << std::endl;
    drawable_lines.push_back(new_drawable_line(B, D, red(0.25f)));
////    drawable_lines.push_back(new_drawable_line(C, E));
////    drawable_lines.push_back(new_drawable_line(A, E));
////    drawable_lines.push_back(new_drawable_line(E, D));

//    std::cout << "Line from origin to A" << std::endl;
//    drawable_lines.push_back(new_drawable_line(origin, A));
//    std::cout << "Line from -x to x" << std::endl;
//    drawable_lines.push_back(new_drawable_line(A, B, cyan()));
//    std::cout << "Line from -y to y" << std::endl;
//    drawable_lines.push_back(new_drawable_line(C, D, magenta()));
//    std::cout << "Line from -z to z" << std::endl;
//    drawable_lines.push_back(new_drawable_line(E, F, yellow()));

    for(auto l : drawable_lines) {
        cubeGeode->addDrawable(l);
    }
//    std::vector<osg::ShapeDrawable*> drawable_planes;
//    drawable_planes.push_back(new_drawable_plane(origin, D, A));
//    for(auto pl : drawable_planes) {
//        cubeGeode->addDrawable(pl);
//    }
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

void motor_orbits()
{
    // The following primitives were inspired by ganja.js
    // https://enkimute.github.io/ganja.js/examples/coffeeshop.html#pga3d_motor_orbits

    // The orbits of parametrised motor variables :
    // these are functions that take one parameter 0<x<1 and return a motor.
    // For a function of 1 parameter, the orbit is rendered as a curve.

    //  var circle  = (BV,r)=>x=>Math.E**(Math.PI*x*BV)*(1+r*.5e01),
    //      segment = (BV)=>x=>1+x*0.5*BV;

    // The product of two such parametrised 1D orbits is a 2D manifold :
    // You can simply multiply the two 1-parameter orbits to
    // make a 2-parameter manifold.

    // this allows us to multiply two circles into a torus or sphere,
    // two segments into a plane, or a segment and a circle to produce
    // a disk, cone or cylinder.

    //  var torus    = (r1,r2)=>circle(1e12,r1)*circle(1e13,r2),
    //      sphere   = (r)=>circle(1e13,0)*circle(1e12,r),
    //      plane    = (x,y)=>segment(x*1e02)*segment(y*1e03),
    //      cylinder = (r,l)=>segment(l*1e02)*circle(1e13,r),
    //      disk     = (r)=>circle(1e12,0)*segment(r*1e01),
    //      cone     = (r,l)=>circle(1e12,0)*segment(r*1e01-l*1e03);
}

/*
scene()
{

    // Visualisation with OpenSceneGraph
    osg::Group* sceneRoot = new osg::Group;

    osg::Box* cube = new osg::Box(osg::Vec3(0.0, 0.0, 0.0), 0.2f);
    osg::ShapeDrawable* cubeDrawable = new osg::ShapeDrawable(cube);
    osg::Geode* cubeGeode = new osg::Geode();
    cubeGeode->addDrawable(cubeDrawable);
    osg::PositionAttitudeTransform* cubeTransform = new osg::PositionAttitudeTransform();
    cubeTransform->addChild(cubeGeode);
    sceneRoot->addChild(cubeTransform);

    osg::Sphere* s1Sphere = new osg::Sphere(osg::Vec3(0.0, 0.0, 0.0), 0.1f);
    osg::ShapeDrawable* s1SphereDrawable = new osg::ShapeDrawable(s1Sphere);
    osg::Geode* s1SphereGeode = new osg::Geode();
    s1SphereGeode->addDrawable(s1SphereDrawable);
    osg::PositionAttitudeTransform* s1SphereTransform = new osg::PositionAttitudeTransform();
    s1SphereTransform->addChild(s1SphereGeode);
    sceneRoot->addChild(s1SphereTransform);

    osg::Sphere* s2Sphere = new osg::Sphere(osg::Vec3(0.0, 0.0, 0.0), 0.1f);
    osg::ShapeDrawable* s2SphereDrawable = new osg::ShapeDrawable(s2Sphere);
    osg::Geode* s2SphereGeode = new osg::Geode();
    s2SphereGeode->addDrawable(s2SphereDrawable);
    osg::PositionAttitudeTransform* s2SphereTransform = new osg::PositionAttitudeTransform();
    s2SphereTransform->addChild(s2SphereGeode);
    sceneRoot->addChild(s2SphereTransform);

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

    // Animation loop: Integration of equations of motion. Generelly Euler-backwards, exception: implicit Euler-forward
    // for solving Euler's equations
    double frameTime = 0.0;
    double sumFrameTime = 0.0;
    double minFrameTime = 0.0;
    double timer = 0.0;
    unsigned int counter = 0;
    while(!viewer.done()) {
        osg::Timer_t startFrameTick = osg::Timer::instance()->tick();

        // Displacement propagation
        // auto dD = part<0x00, 0x03, 0x05, 0x06, 0x09, 0x0a, 0x0c, 0x0f, 0x11, 0x12, 0x14, 0x17>(D*V_b*0.5);
        auto dD = part_type<D_type>(D * V_b * 0.5);
        D = D + dD * frameTime;

        // Generalized trigonometric formula for the tangent of a half angle (after Hestenes)
        auto B_s1 = part_type<S_type>(
            grade<2>(grade<2>((~D) * D_s1) * (!(one + grade<0>((~D) * D_s1) + grade<4>((~D) * D_s1)))));
        auto B_s2 = part_type<S_type>(
            grade<2>(grade<2>((~D) * D_s2) * (!(one + grade<0>((~D) * D_s2) + grade<4>((~D) * D_s2)))));

        // Force law of spring (not necessarily linear)
        auto F_s1_b = B_s1 * k_s1;
        auto F_s2_b = B_s2 * k_s2;
        // Force law of damping
        auto F_d = V_b * k_d * (-1.0);
        // Gravity acting on body
        auto F_g = part_type<S_type>((!part<0, 3, 5, 6>(D)) * einf * (e3 * (-9.81)) * part<0, 3, 5, 6>(D));
        // Resultant force wrench
        auto F_b = eval(F_s1_b + F_s2_b + F_d + F_g);

        //--- Start: velocity propagation due to force laws ---
        // Torque part of force wrench
        auto t_b = eval(i * part<0x03, 0x05, 0x06>(F_b));
        // Linear force part of force wrench
        auto f_b = grade<1>(part<0x09, 0x0a, 0x0c, 0x11, 0x12, 0x14>(F_b) * e0);

        // Linear velocity part of velocity twist propagation
        auto v = grade<1>(part<0x09, 0x0a, 0x0c, 0x11, 0x12, 0x14>(V_b) * e0) + f_b * (frameTime / M);

        // Angular velocity part of velocity twist propagation (implicit Euler-forward, solving Euler's equations)
        auto oldOm = eval((~i) * (-1.0) * part<0x03, 0x05, 0x06>(V_b));
        auto om = oldOm;
        decltype(om) prevOm;
        double maxError = 1e-5;
        // do {
        for(int j = 0; j < 5; ++j) {
            prevOm = om;
            om[0] = oldOm[0] + (t_b[0] - (In_3 - In_2) * om[1] * om[2]) / In_1 * frameTime;
            om[1] = oldOm[1] + (t_b[1] - (In_1 - In_3) * om[2] * om[0]) / In_2 * frameTime;
            om[2] = oldOm[2] + (t_b[2] - (In_2 - In_1) * om[0] * om[1]) / In_3 * frameTime;
        }
        //} while((magnitude(om-prevOm).element<0x00>()) > maxError);

        // Combining velocity propagations to velocity twist
        V_b = i * om * (-1.0) + einf * v;
        //--- End: velocity propagation due to force laws ---

        // Updating new position of cube
        auto p_m = eval(grade<1>(D * e0 * (~D)));
        cubeTransform->setPosition(osg::Vec3(p_m[0], p_m[1], p_m[2]));
        cubeTransform->setAttitude(osg::Quat(-D[3], D[2], -D[1], D[0]));

        // Updating actuating points
        auto T_s1 = (one + einf * (e1 * 1.0 + e2 * 0.0 + e3 * 0.5 * (1.0 + sin(5.0 * timer))) * 0.5);
        D_s1 = T_s1 * R_s1;
        auto p_s1 = eval(grade<1>(D_s1 * e0 * (~D_s1)));
        s1SphereTransform->setPosition(osg::Vec3(p_s1[0], p_s1[1], p_s1[2]));

        auto T_s2 = (one + einf * (e1 * 0.0 + e2 * sin(3.0 * timer) + e3 * 0.0) * 0.5);
        D_s2 = T_s2 * R_s2;
        auto p_s2 = eval(grade<1>(D_s2 * e0 * (~D_s2)));
        s2SphereTransform->setPosition(osg::Vec3(p_s2[0], p_s2[1], p_s2[2]));

        viewer.frame();

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

        if(frameTime < minFrameTime)
            OpenThreads::Thread::microSleep(static_cast<unsigned int>(1000000.0 * (minFrameTime - frameTime)));
    }
}

int main()
{

    auto a_line = pga3::line(2, 3, 4, 5, 6, 7);
    auto a_plane = pga3::plane(2, 3, 4, 7);
    auto a_point = pga3::point(2, 3, 4);
    auto an_ideal_point = pga3::ideal_point(2, 3, 4);
    auto another_point = pga3::point(4, 7, 8) * 2.;
    auto A = pga3::point(0, -1, 0), B = pga3::point(1, 1, -1), C = pga3::point(-1, 1, -1), D = pga3::point(1, 1, 1),
         E = pga3::point(-1, 1, 1);
    auto centroid = A + B + C + D + E;

    std::cout << "centroid: " << centroid << " normalized centroid: " << normalize(centroid) << std::endl;
    std::cout << "inverse centroid: " << ::inverse(normalize(centroid)) << std::endl;

    std::cout << "A line from points(centroid,C+E): " << pga3::line_from_points(centroid, C + E) << std::endl;
    std::cout << "Sum of points: " << normalize(a_point + another_point) << std::endl;
    std::cout << "A point doubled then normalized: " << normalize(2. * a_point) << std::endl;
    std::cout << "A point normalized: " << normalize(a_point) << std::endl;
    std::cout << "Another point normalized: " << normalize(another_point) << std::endl;
    std::cout << "An ideal point dualized: " << ::dual(an_ideal_point) << std::endl;
    // Next line fails to compile since can't normalize an ideal point
    //    std::cout << "An ideal point normalized: " << normalize(an_ideal_point) << std::endl;
    std::cout << "A point dualized: " << ::dual(a_point) << std::endl;
    std::cout << "A line normalized: " << normalize(a_line) << std::endl;
    std::cout << "A line dualized: " << ::dual(a_line) << std::endl;
    std::cout << "A plane normalized: " << normalize(a_plane) << std::endl;
    std::cout << "A plane dualized, then normalized: " << normalize(::dual(a_plane)) << std::endl;
    std::cout << "A line from points: " << pga3::line_from_points(a_point, another_point) << std::endl;
    std::cout << "Squared line from points: "
              << pga3::line_from_points(a_point, another_point) * pga3::line_from_points(a_point, another_point)
              << std::endl;
    std::cout << "Squared line from points: "
              << pga3::line_from_points(a_point, another_point) * pga3::line_from_points(a_point, another_point)
              << std::endl;
    std::cout << "A plane from points: " << pga3::plane_from_points(A, B, C) << std::endl;
    std::cout << "A line from planes: "
              << pga3::line_from_planes(pga3::plane_from_points(A, B, C), pga3::plane_from_points(B, C, D))
              << std::endl;
    std::cout << "A point from line and plane: "
              << pga3::point_from_line_and_plane(
                     pga3::line_from_points(centroid, C + E), pga3::plane_from_points(A, D, B))
              << std::endl;
    std::cout << "sum of lines:" << (pga3::line_from_points(B, D) + pga3::line_from_points(C, E)) << std::endl;
    std::cout << "motor sandwich applied to a line:"
              << pga3::sandwich(
                     pga3::line_from_points(C, E), pga3::motor(pga3::line_from_points(centroid, C + E), 2, M_PI / 3.))
              << std::endl;
    std::cout << "motor2 sandwich applied to a line:"
              << pga3::sandwich(
                     pga3::line_from_points(C, E), pga3::screw2(pga3::line_from_points(centroid, C + E), 2, M_PI / 3.))
              << std::endl;
    std::cout << "motor sandwich applied to point:" << A << " "
              << (pga3::sandwich(A, pga3::motor(pga3::line_from_points(centroid, C + E), 2, M_PI / 3.))) << std::endl;
    std::cout << "motor2 sandwich applied to point:" << A << " "
              << (pga3::sandwich(A, pga3::screw2(pga3::line_from_points(centroid, C + E), 2, M_PI / 3.))) << std::endl;
    std::cout << "translator sandwich applied to point:" << A << " "
              << (pga3::sandwich(A, pga3::translator(pga3::line_from_points(centroid, C + E), 1))) << std::endl;
    std::cout << "vee applied to line:" << pga3::vee(pga3::E0, 0.5 * normalize(pga3::line_from_points(centroid, C + E)))
              << std::endl;

    //    std::cout << "Pseudoscalar normalized: " << normalize(-3*pga3::I) << std::endl;
}
*/