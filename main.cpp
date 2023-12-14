#include <iostream>
#include <cmath>
#include <ostream>

#include "gaalet.h"
#include "pga3.h"
#include "pga3_point.h"
#include "pga3_line.h"
#include "pga3_ops.h"
#include "pga3_normalize.h"
#include "pga3_logarithm.h"
#include "pga3_exponential.h"
#include "Eigen/Core"

using namespace gaalet;
using namespace Eigen;

void rotation_around_x_axis(double angle, Matrix4d& rotation_)
{
    // Return the 4x4 rotation matrix for a rotation around the x-axis

    double c = cos(angle);
    double s = sin(angle);
		Matrix2d rotation;
    rotation_(0, 0) = 1.;
    rotation_(0, 1) = 0.;
    rotation_(0, 2) = 0.;
    rotation_(0, 3) = 0.;
    rotation_(1, 0) = 0.;
    rotation_(1, 1) = c;
    rotation_(1, 2) = -s;
    rotation_(1, 3) = 0.;
    rotation_(2, 0) = 0.;
    rotation_(2, 1) = s;
    rotation_(2, 2) = c;
    rotation_(2, 3) = 0.;
    rotation_(3, 0) = 0.;
    rotation_(3, 1) = 0.;
    rotation_(3, 2) = 0.;
    rotation_(3, 3) = 1.;
}

void rotation_around_y_axis(double angle, Matrix4d& rotation_)
{
// Return the 4x4 rotation matrix for a rotation around the y-axis

    double c = cos(angle);
    double s = sin(angle);
    rotation_(0, 0) = c;
    rotation_(0, 1) = 0.;
    rotation_(0, 2) = s;
    rotation_(0, 3) = 0.;
    rotation_(1, 0) = 0;
    rotation_(1, 1) = 1.;
    rotation_(1, 2) = 0.;
    rotation_(1, 3) = 0.;
    rotation_(2, 0) = -s;
    rotation_(2, 1) = 0.;
    rotation_(2, 2) = c;
    rotation_(2, 3) = 0.;
    rotation_(3, 0) = 0.;
    rotation_(3, 1) = 0.;
    rotation_(3, 2) = 0.;
    rotation_(3, 3) = 1.;
}

void rotation_around_z_axis(double angle, Matrix4d& rotation_)
{
// Return the 4x4 rotation matrix for a rotation around the z-axis

    double c = cos(angle);
    double s = sin(angle);
    rotation_(0, 0) = c;
    rotation_(0, 1) = -s;
    rotation_(0, 2) = 0.;
    rotation_(0, 3) = 0.;
    rotation_(1, 0) = s;
    rotation_(1, 1) = c;
    rotation_(1, 2) = 0.;
    rotation_(1, 3) = 0.;
    rotation_(2, 0) = 0.;
    rotation_(2, 1) = 0.;
    rotation_(2, 2) = 1.;
    rotation_(2, 3) = 0.;
    rotation_(3, 0) = 0.;
    rotation_(3, 1) = 0.;
    rotation_(3, 2) = 0.;
    rotation_(3, 3) = 1.;
}

void translation(Vector3d t, Matrix4d& t_matrix_)
{
    // Return the 4x4 translation matrix
		t_matrix_.setIdentity();
    t_matrix_(0, 3) = t[0];
    t_matrix_(1, 3) = t[1];
    t_matrix_(2, 3) = t[2];
}

double degrees2radians(double x) {
    return x/180.0 * M_PI;
}

int main() {


    double x_rotation_angle = 15.0;
    double y_rotation_angle = 32.0;
    double z_rotation_angle = 68.0;

		Vector4d test_point;
    test_point(0) = 1.0;
    test_point(1) = 2.0;
    test_point(2) = 3.0;
    test_point(3) = 1.0;

    pga3::Point_t pga_test_pt = pga3::make_point(test_point(0), test_point(1), test_point(2));

    pga3::Point_t pga_test_pt2 = pga3::make_point(3., 4.5, -16);
    std::cout << "pga_test_pt2: " << pga3::Point(pga_test_pt2) << std::endl;

		Matrix4d rm_x;
    rotation_around_x_axis(degrees2radians(x_rotation_angle), rm_x);
    std::cout << "rm_x: " << std::endl << rm_x * test_point << std::endl;
		Matrix4d rm_y;
    rotation_around_y_axis(degrees2radians(y_rotation_angle), rm_y);
    std::cout << "rm_y: " << std::endl << rm_y * test_point << std::endl;
		Matrix4d rm_z;
    rotation_around_z_axis(degrees2radians(z_rotation_angle), rm_z);
    std::cout << "rm_z * (1, 1, 1).T: " << std::endl << rm_z * test_point << std::endl;

		Matrix4d trans;
    translation(Vector3d (3., 7., 9.), trans);
    std::cout << trans << std::endl;

		Matrix4d rm = rm_x * rm_y * rm_z;
    std::cout << "rotated with rm: " << std::endl << rm * test_point << std::endl;

		Matrix4d tr_rm = trans * rm_x * rm_y * rm_z;
    std::cout << "translated and rotated with rm: " << std::endl << tr_rm * test_point << std::endl;

    auto r_x = pga3::rotor(pga3::k, degrees2radians(x_rotation_angle));
    std::cout << "r_x: " << eval(r_x) << std::endl;
    std::cout << "Point(1,1,1) rotated around x-axis by " << x_rotation_angle << " degrees: "
              << pga3::Point(r_x * pga_test_pt * (~r_x)) << std::endl;
    auto r_y = pga3::rotor(pga3::j, degrees2radians(y_rotation_angle));
    std::cout << "Point(1,1,1) rotated around y-axis by " << y_rotation_angle << " degrees: "
              << pga3::Point(r_y * pga_test_pt * (~r_y)) << std::endl;
    auto r_z = pga3::rotor(pga3::i, degrees2radians(z_rotation_angle));
    std::cout << "Point(1,1,1) rotated around z-axis by " << z_rotation_angle << " degrees: "
              << pga3::Point(r_z * pga_test_pt * (~r_z)) << std::endl;

    auto t_r = pga3::translator(pga3::line_from_points(pga3::E0, pga3::make_point(3., 7., 9.)),
            sqrt(3.*3.+7.*7.+9.*9.));
    std::cout << "pga translator: " << eval(t_r) << std::endl;

    auto r_xyz = r_x * r_y * r_z;
    auto t_r_xyz = t_r * r_x * r_y * r_z;
    auto rp = pga3::Point(pga3::sandwich(pga_test_pt, r_xyz));
    auto tr_rp = pga3::Point(pga3::sandwich(pga_test_pt, t_r_xyz));

    std::cout << "Point(1,1,1) rotated around rotations around z, then y, then x axes "
              << rp << std::endl;

    std::cout << "Point(1,1,1) rotated around rotations around z, then y, then x axes then translated "
              << tr_rp << std::endl;

//    std::cout << "exp(5.): " << eval(pga3::exp(5.0 * pga3::one)) << std::endl;

    auto log_r_xyz = pga3::log(r_xyz);
    auto log_r_x = pga3::log(r_x);
    auto log_r_y = pga3::log(r_y);
    auto log_r_z = pga3::log(r_z);

    std::cout << "Composite rotor" << std::endl;
    std::cout << eval(r_xyz) << std::endl << std::endl;
    std::cout << "Logarithms" << std::endl;
    std::cout << "log_r_x: " << eval(log_r_x) << "  exp(log_r_x): " << eval(pga3::exp(log_r_x)) << "  r_x: " << eval(r_x) << std::endl;
    std::cout << "log_r_y: " << eval(log_r_y) << "  exp(log_r_y): " << eval(pga3::exp(log_r_y)) << "  r_y: " << eval(r_y)  << std::endl;
    std::cout << "log_r_z: " << eval(log_r_z) << "  exp(log_r_z): " << eval(pga3::exp(log_r_z)) << "  r_z: " << eval(r_z)  << std::endl;
    std::cout << "log_r_xyz: " << eval(log_r_xyz) << std::endl << "exp(log_r_xyz): " << eval(pga3::exp(log_r_xyz)) << std::endl;
    std::cout << "exp(log_r_x) * exp(log_r_y) * exp(log_r_z): "
              <<   eval(pga3::exp(log_r_x) * pga3::exp(log_r_y) * pga3::exp(log_r_z)) << std::endl;
    auto RXYZ = eval(pga3::exp(log_r_x) * pga3::exp(log_r_y) * pga3::exp(log_r_z));
    std::cout << "RXYZ = exp(log_r_x) * exp(log_r_y) * exp(log_r_z): "
              <<  RXYZ  << std::endl;
    std::cout << "(RXYZ * (~RXYZ)): " << (RXYZ * (~RXYZ)) << std::endl;
    std::cout << "r_x * r_y : " <<  eval(r_x * r_y)  << std::endl;
    std::cout << "log(r_x * r_y) : " <<  eval(pga3::log(r_x * r_y))  << std::endl;
    std::cout << "exp(log(r_x) + log(r_y)) : " <<  eval(pga3::exp(pga3::log(r_x) + pga3::log(r_y)))  << std::endl;
    std::cout << "study_number(r_x * r_y) : " <<  pga3::bivector_axis(r_x * r_y).first.first << " + "
              <<  pga3::bivector_axis(r_x * r_y).first.second << std::endl;
    std::cout << "bivector_axis(r_x * r_y) : " <<  eval(pga3::bivector_axis(r_x * r_y).second)  << std::endl;
    std::cout << "study_number(r_x * r_y * r_z) : " <<  pga3::bivector_axis(r_x * r_y * r_z).first.first << " + "
              <<  pga3::bivector_axis(r_x * r_y * r_z).first.second << std::endl;
    std::cout << "bivector_axis(r_x * r_y * r_z) : " <<  eval(pga3::bivector_axis(r_x * r_y * r_z).second)  << std::endl;
    std::cout << "log(bivector_axis(r_x * r_y * r_z)) : " <<  eval(pga3::log(pga3::bivector_axis(r_x * r_y * r_z).second))  << std::endl;
    std::cout << "sn_a * bivector_axis(r_x * r_y * r_z) : " <<  eval(pga3::bivector_axis(r_x * r_y * r_z).first.first * pga3::bivector_axis(r_x * r_y * r_z).second)  << std::endl;
    std::cout << "study_number(log(r_x * r_y * r_z)) : " <<  pga3::bivector_axis(pga3::log(r_x * r_y * r_z)).first.first << " + "
            <<  pga3::bivector_axis(pga3::log(r_x * r_y * r_z)).first.second << std::endl;
    std::cout << "bivector_axis(log(r_x * r_y * r_z)) : " <<  eval(pga3::bivector_axis(pga3::log(r_x * r_y * r_z)).second)  << std::endl;
    std::cout << "sn_a * bivector_axis(log(r_x * r_y * r_z)) : " <<  eval(pga3::bivector_axis(pga3::log(r_x * r_y * r_z)).first.first * pga3::bivector_axis(pga3::log(r_x * r_y * r_z)).second)  << std::endl;

    auto a_line = pga3::line_from_points(pga3::make_point(1., 1., 1.),
            pga3::make_point(2., 3., 4.));
    std::cout << "Line from (1,1,1) to (2, 3, 4) " << eval(a_line) << std::endl;

    auto l = pga3::Line(a_line);
    std::cout << "Line from (1,1,1) to (2, 3, 4) " << l << ", " << std::endl;
    std::cout << "Line from (1,1,1) to (2, 3, 4) " << eval(pga3::normalize(a_line)) << ", " << std::endl;
    std::cout << "Line from (1,1,1) to (2, 3, 4) " << l.normalized() << std::endl;

    auto rotated_line = pga3::sandwich(a_line, r_xyz);
    std::cout << "Line from (1,1,1) to (2, 3, 4) rotated around around x, then y, then z axes "
              << eval(::grade<2>(rotated_line)) << std::endl;

}