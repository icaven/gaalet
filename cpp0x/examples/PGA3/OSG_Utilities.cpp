// Utility functions to help draw scenes using OpenSceneGraph

#include "OSG_Utilities.h"


// Functions to create drawable primitives

/// Draw a sphere for a point
osg::ShapeDrawable* new_drawable_point(const pga3::Point_t& p, const osg::Vec4& colour,
                                       float point_radius)
{
    osg::Sphere* sphere = new osg::Sphere(Vec3(p), point_radius);
    osg::ShapeDrawable* drawable = new osg::ShapeDrawable(sphere);
    drawable->setColor(colour);
    return drawable;
}


/// Draw a cylinder for a line
osg::ShapeDrawable*
new_drawable_line(const pga3::Point_t& start_pt, const pga3::Point_t& end_pt, 
                  const osg::Vec4& colour,
                  const float line_thickness)
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
        cylinder->setCenter(cylinder->getCenter() + osg::Vec3(0., 0., -length*0.5f));
    }
    else if (isclose(angle, M_PI)) {
        // Parallel to the Z-axis, just shift along that axis in the negative direction
        cylinder->setCenter(cylinder->getCenter() + osg::Vec3(0., 0., length*0.5f));
    }
    else {
        // Determine the plane that the Z-axis and the line forms and then compute the perpendicular
        // to that plane; the angle will correspond to the angle around that perpendicular
        auto end_of_cylinder = normalize(pga3::sandwich(start_pt, pga3::translator(pga3::k, length)));
        auto perpendicular_to_plane = pga3::normal_to_plane(start_pt,  end_pt, end_of_cylinder);

        // Rotate the cylinder to be in the direction of the line
        osg::Quat q = Quat(M_PI/2.0+angle/2.0, perpendicular_to_plane);
        cylinder->setRotation(q);
    
        // Translate the center point
        auto t = pga3::translator(line, length*0.5);
        auto new_center = pga3::sandwich(start_pt, t);
        cylinder->setCenter(Vec3(new_center));
    }

    osg::ShapeDrawable* drawable = new osg::ShapeDrawable(cylinder);
    drawable->setColor(colour);
    return drawable;
}


/// Arrow shape
osg::CompositeShape*
new_arrow(const pga3::Point_t& start_pt, const pga3::Point_t& end_pt, 
          const float line_thickness)
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
        shaft->setCenter(shaft->getCenter() + osg::Vec3(0., 0., -shaft_length*0.5f));
        arrow_head->setRotation(osg::Quat(M_PI, osg::Vec3d(1., 0., 0.))); // Reflect the arrow head
        arrow_head->setCenter(shaft->getCenter() + 
            osg::Vec3(0., 0., fmin(-(length-head_offset)*0.5f, -arrow_head_length)));
    }
    else if (isclose(angle, M_PI)) {
        // Parallel to the Z-axis, just shift along that axis in the negative direction
        shaft->setCenter(shaft->getCenter() + osg::Vec3(0., 0., shaft_length*0.5f));
        arrow_head->setCenter(shaft->getCenter() + 
            osg::Vec3(0., 0., fmax((length-head_offset)*0.5f, arrow_head_length)));
   }
    else {
        // Determine the plane that the Z-axis and the line forms and then compute the perpendicular
        // to that plane; the angle will correspond to the angle around that perpendicular
        auto end_of_shaft = normalize(pga3::sandwich(start_pt, pga3::translator(pga3::k, shaft_length)));
        auto perpendicular_to_plane = pga3::normal_to_plane(start_pt, end_pt, end_of_shaft);

        // Rotate the shaft to be in the direction of the line
        osg::Quat q = Quat(M_PI/2.+angle/2.0, perpendicular_to_plane);
        shaft->setRotation(q);

        // Translate the center point of the arrow shaft
        auto t = pga3::translator(line, shaft_length*0.5);
        auto new_center = pga3::sandwich(start_pt, t);
        shaft->setCenter(Vec3(new_center));
        
        // Rotate and translate the arrow head
        arrow_head->setRotation(q);
//        auto arrow_translator = pga3::translator(line, (double) fmax(length-head_offset, arrow_head_length));
        auto arrow_translator = pga3::translator(line, (double) fmax(length-head_offset, arrow_head_length-head_offset));
        auto arrow_center = pga3::sandwich(start_pt, arrow_translator);
        arrow_head->setCenter(Vec3(arrow_center));
    }

    osg::CompositeShape* arrow = new osg::CompositeShape();
    arrow->addChild(shaft);
    arrow->addChild(arrow_head);
    
    return arrow;
}


/// Drawable arrow
osg::ShapeDrawable*
new_drawable_arrow(const pga3::Point_t& start_pt, const pga3::Point_t& end_pt, 
                  const osg::Vec4& colour,
                  const float line_thickness)
{
    osg::CompositeShape* arrow = new_arrow(start_pt, end_pt, line_thickness);
    osg::ShapeDrawable* drawable = new osg::ShapeDrawable(arrow);
    drawable->setColor(colour);
    return drawable;
}


/// Drawable arrow
osg::ShapeDrawable*
new_drawable_arrow(const pga3::Point_t& start_pt, const pga3::Line_t& line, 
                  const osg::Vec4& colour,
                  const float line_thickness)
{
    osg::CompositeShape* arrow = new_arrow(start_pt, start_pt + pga3::make_ideal_point(line), line_thickness);
    osg::ShapeDrawable* drawable = new osg::ShapeDrawable(arrow);
    drawable->setColor(colour);
    return drawable;
}

osg::ShapeDrawable*
new_drawable_plane(const pga3::Point_t& p1, const pga3::Point_t& p2, const pga3::Point_t& p3,
                   const osg::Vec4& colour, const float plane_thickness)
{
    // Determine the plane, its normal vector, and the point at the end of the normal
    auto perpendicular_to_plane = pga3::normal_to_plane(p1, p2, p3);
    auto translate_along_normal = pga3::translator(perpendicular_to_plane, 1.0);
    auto end_of_normal = pga3::sandwich(p1, translate_along_normal);
    
    // The box that will be used to represent the plane is initially oriented so that
    // its thickness is in the Y-axis.  There a rotation will need to be computed to
    // rotate the box into place.  The axis of this rotation is the perpendicular to the
    // plane formed by the p1 point (the initial center of the box), the end of the normal
    // to the original desired plane, and the point that is unit distance "upwards" from
    // the p1 point.  So there are two planes involved: the desired plane represented by a
    // box, and the plane of rotation of that box.
    auto up = pga3::translator(pga3::j, 1.0);
    auto p1_upwards = pga3::sandwich(p1, up);

    auto perpendicular_to_plane_with_normal = pga3::normal_to_plane(p1,  end_of_normal, p1_upwards);
    auto translate_along_normal_rotation = pga3::translator(perpendicular_to_plane_with_normal, 1.0);
    auto end_of_normal_to_rotation_plane = pga3::sandwich(p1, translate_along_normal_rotation);

    auto bisecting_line = pga3::line_from_points(p1, normalize(p2 + p3));
    double x_length = eval(::magnitude(bisecting_line));
    pga3::Line_t other_line = pga3::line_from_points(p2, p3);
    double z_length = eval(::magnitude(other_line));
    auto towards_bisecting_point = pga3::translator(bisecting_line, x_length/2.0) * pga3::translator(other_line, z_length/4.0);
    auto plane_center = pga3::sandwich(p1, towards_bisecting_point);

    double angle = acos((pga3::j & perpendicular_to_plane).template element<0>());

    osg::Box* box_as_plane = new osg::Box(Vec3(p1),float(x_length), plane_thickness, float(z_length));
    osg::Quat q = Quat(angle/2.0, perpendicular_to_plane_with_normal);
    box_as_plane->setRotation(q);
    box_as_plane->setCenter(Vec3(plane_center));
    osg::ShapeDrawable* drawable = new osg::ShapeDrawable(box_as_plane);

    // For debugging, it is helpful to draw the normals together with the rotated box,
    // instead of just the box
//    osg::CompositeShape* plane_and_normal= new osg::CompositeShape();
//    plane_and_normal->addChild(box_as_plane);
//    plane_and_normal->addChild(new osg::Sphere(Vec3(plane_center), DEFAULT_RADIUS_OF_DRAWN_POINT));
//    plane_and_normal->addChild(new osg::Sphere(Vec3(normalize(p2 + p3)), DEFAULT_RADIUS_OF_DRAWN_POINT));
//    plane_and_normal->addChild(new_arrow(p1, normalize(p2 + p3), DEFAULT_RADIUS_OF_DRAWN_POINT/3.));
//
//    plane_and_normal->addChild(new_arrow(p1, end_of_normal));
//    plane_and_normal->addChild(new_arrow(p1, end_of_normal_to_rotation_plane));
//
//    osg::ShapeDrawable* drawable = new osg::ShapeDrawable(plane_and_normal);
    
    
    drawable->setColor(colour);
    return drawable;
}

osg::ShapeDrawable*
new_drawable_triangle(const pga3::Point_t& p1, const pga3::Point_t& p2, const pga3::Point_t& p3,
                      const osg::Vec4& colour, bool draw_normal)
{
    osg::TriangleMesh* triangle_as_plane = new osg::TriangleMesh();
    osg::Vec3Array* vertices = new osg::Vec3Array();
    vertices->push_back(Vec3(p1));
    vertices->push_back(Vec3(p2));
    vertices->push_back(Vec3(p3));
    triangle_as_plane->setVertices(vertices);
    std::vector<int> triangle_indices;
    triangle_indices.push_back(0);
    triangle_indices.push_back(1);
    triangle_indices.push_back(2);
    osg::IntArray* indices = new osg::IntArray(triangle_indices.begin(), triangle_indices.end());
    triangle_as_plane->setIndices(indices);

    osg::ShapeDrawable* drawable;
    if (draw_normal) {
        // For debugging, it is helpful to draw the normal together with the triangle
        auto translate_along_normal = pga3::translator(pga3::normal_to_plane(p1, p2, p3), 1.0);
        auto triangle_centroid = eval(::normalize(p1 + p2 + p3));
        auto end_of_normal = pga3::sandwich(triangle_centroid, translate_along_normal);

        osg::CompositeShape* triangle_and_normal= new osg::CompositeShape();
        triangle_and_normal->addChild(triangle_as_plane);
        triangle_and_normal->addChild(new_arrow(triangle_centroid, end_of_normal));
        drawable = new osg::ShapeDrawable(triangle_and_normal);
    }
    else {
        drawable = new osg::ShapeDrawable(triangle_as_plane);
    }

    drawable->setColor(colour);
    return drawable;
    
}