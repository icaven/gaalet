
#include <cmath>
#include <vector>

#include "pga3_streaming.h"

#include "gaalet.h"
#include "pga3.h"
#include "pga3_ops.h"
#include "pga3_norm.h"
#include "pga3_normalize.h"
#include "pga3_dual.h"
#include "pga3_point.h"
#include "pga3_line.h"
#include "pga3_plane.h"
#include "pga3_utility.h"

using namespace gaalet;

template <class MV>
void print_basis_info(const MV& p, const std::string& name) {
    std::cout << name << ": " << p 
            << " polar of " << name << ": " << eval(pga3::polar(p))
            << " dual of " << name << ": " << eval(pga3::dual(p))
            << " " << name << " ^ dual(" << name << ") " << eval(p ^ pga3::dual(p))
            << " square of " << name << ": " << eval(p * p) // + 0. * pga3::one)
            << " normalized " << name << ": " << eval(pga3::normalize(p)) // + 0. * pga3::one)
            << std::endl;
}


void print_point_info(const pga3::Point_t& p, const std::string& name) {
    std::cout << name << ": " << eval(p)
              << " dual of " << name << ": " << eval(pga3::dual(p)) << std::endl;
    std::cout << "normalized " << name << ": " << eval(pga3::normalize(p)) << std::endl;
}


void print_line_info(const pga3::Line_t& l, const std::string& name) {
    std::cout << name << ": " << l 
              << " polar of " << name << ": " << eval(pga3::polar(l))
              << " dual of " << name << ": " << eval(pga3::dual(l)) << std::endl;
    std::cout << "normalized " << name << ": " << eval(pga3::normalize(l)) << std::endl;
}

void print_cayley_table()
{
    typedef  pga3::space::mv<0,
            pga3::vector_conf(0), pga3::vector_conf(1), pga3::vector_conf(2), pga3::vector_conf(3),
            blade_conf(0, 1), blade_conf(0, 2), blade_conf(0, 3), blade_conf(1, 2), blade_conf(3, 1), blade_conf(2, 3),
            blade_conf(0, 2, 1), blade_conf(0, 1, 3), blade_conf(0, 3, 2), blade_conf(1, 2, 3),
            pga3::pseudoscalar_conf>::type all_conf_list;

    std::vector<all_conf_list> basis_elements = {pga3::one, pga3::e0, pga3::e1, pga3::e2, pga3::e3,
                                                 pga3::e01, pga3::e02, pga3::e03, pga3::e12, pga3::e31,
                                                 pga3::e23, pga3::e021, pga3::e013, pga3::e032, pga3::e123, pga3::I};

    std::cout << "The basis elements and their duals:" << std::endl;
    for (auto element1: basis_elements) {
        std::cout << "\"" << std::setw(8) << std::left << element1 << "\" "
                  << "dual: \"" << std::setw(8) << std::left << eval(pga3::dual(element1)) << "\" " << std::endl;
    }
    std::cout << std::endl;

    std::cout << "The Cayley table for PGA3:" << std::endl;
    for (auto element1: basis_elements) {
        for (auto element2: basis_elements) {
            std::cout << "\"" << eval(element1 * element2) << "\" ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

}

void print_info_all_basis()
{
    std::cout << "total dimensions: " << pga3::space::algebra::metric::dimension << std::endl;
    std::cout << "signature_bitmap: " << std::setbase(16) << pga3::space::algebra::metric::signature_bitmap << std::endl;
    std::cout << "euclidean_bitmap: " << std::setbase(16) << pga3::space::algebra::metric::euclidean_bitmap << std::endl;
    std::cout << "degenerate_bitmap: " << std::setbase(16) << pga3::space::algebra::metric::degenerate_bitmap << std::endl;

    print_basis_info(pga3::one, "1");
    print_basis_info(pga3::e0, "e0");
    print_basis_info(pga3::e1, "e1");
    print_basis_info(pga3::e2, "e2");
    print_basis_info(pga3::e3, "e3");
    print_basis_info(pga3::e01, "e01");
    print_basis_info(pga3::e02, "e02");
    print_basis_info(pga3::e03, "e03");
    print_basis_info(pga3::e12, "e12");
    print_basis_info(pga3::e31, "e31");
    print_basis_info(pga3::e23, "e23");
    print_basis_info(pga3::e021, "e021");
    print_basis_info(pga3::e013, "e013");
    print_basis_info(pga3::e032, "e032");
    print_basis_info(pga3::e123, "e123");
    print_basis_info(pga3::I, "I");

    // Different names for some of the above entities
    print_basis_info(pga3::i, "i");  
    print_basis_info(pga3::j, "j");
    print_basis_info(pga3::k, "k");
    print_basis_info(pga3::dual_i, "dual_i");
    print_basis_info(pga3::dual_j, "dual_j");
    print_basis_info(pga3::dual_k, "dual_k");
    print_basis_info(pga3::E0, "E0");
    print_basis_info(pga3::E1, "E1");
    print_basis_info(pga3::E2, "E2");
    print_basis_info(pga3::E3, "E3");
}

void print_bivector_info() {
    std::cout << "(pga3::i*pga3::j) " << eval(pga3::i*pga3::j) << std::endl;
    std::cout << "(pga3::i*pga3::k) " << eval(pga3::i*pga3::k) << std::endl;

    std::cout << "(pga3::j*pga3::i) " << eval(pga3::j*pga3::i) << std::endl;
    std::cout << "(pga3::j*pga3::k) " << eval(pga3::j*pga3::k) << std::endl;
    
    std::cout << "(pga3::k*pga3::i) " << eval(pga3::k*pga3::i) << std::endl;
    std::cout << "(pga3::k*pga3::j) " << eval(pga3::k*pga3::j) << std::endl;

}

void print_plane_info() {
    std::cout << eval(pga3::e0 ^ pga3::dual(pga3::e0)) << std::endl;
    std::cout << eval(pga3::e1 ^ pga3::dual(pga3::e1)) << std::endl;
    std::cout << eval(pga3::e2 ^ pga3::dual(pga3::e2)) << std::endl;
    std::cout << eval(pga3::e3 ^ pga3::dual(pga3::e3)) << std::endl;
}

double degrees2radians(double x) {
	return x/180.0 * M_PI;
}

void test_with_points()
{
    std::cout << " " << pga3::E0 << std::endl;

    auto origin = pga3::make_point(0, 0, 0);
    auto A = pga3::make_point(0, -1, 0);
    auto B = pga3::make_point(1, 1, -1);
    auto C = pga3::make_point(-1, 1, -1);
    auto D = pga3::make_point(1, 1, 1);
    auto E = pga3::make_point(-1, 1, 1);

    print_point_info(A, "A");
    print_point_info(B, "B");

    double x_rotation_angle = 15.0;
    double y_rotation_angle = 32.0;
    double z_rotation_angle = 68.0;
    auto test_point = pga3::make_point(1., 1., 1.);
    auto r_x = pga3::rotor(pga3::e2 ^ pga3::e3, degrees2radians(x_rotation_angle));
    std::cout << "Point(1,1,1) rotated around x-axis by " << x_rotation_angle << " degrees: "
              << pga3::Point(r_x * test_point * (~r_x)) << std::endl;
    auto r_y = pga3::rotor(pga3::e1 ^ pga3::e3, degrees2radians(y_rotation_angle));
    std::cout << "Point(1,1,1) rotated around y-axis by " << y_rotation_angle << " degrees: "
              << pga3::Point(r_y * test_point * (~r_y)) << std::endl;
    auto r_z = pga3::rotor(pga3::e1 ^ pga3::e2, degrees2radians(z_rotation_angle));
    std::cout << "Point(1,1,1) rotated around z-axis by " << z_rotation_angle << " degrees: "
              << pga3::Point(r_z * test_point * (~r_z)) << std::endl;
    auto r_zyx = r_z * (r_y * r_x);
    std::cout << "Point(1,1,1) rotated around rotations around z, then y, then x axes "
              << pga3::Point(r_zyx * test_point * (~r_zyx)) << std::endl;
    auto r_xyz = r_x * (r_y * r_z);
    std::cout << "Point(1,1,1) rotated around rotations around x, then y, then z axes "
              << pga3::Point(r_xyz * test_point * (~r_xyz)) << std::endl;

    double x_distance = 2.0;
    auto t_x = pga3::translator(pga3::e23, x_distance);
    auto translated_plane = pga3::sandwich(pga3::e1, t_x);
    std::cout << "Plane(x=0) " << eval(pga3::e1) << " translated on the x-axis "
              << eval(translated_plane) << "then meet with point on that plane"
              << eval(translated_plane ^ pga3::make_point(x_distance, 3, 4)) << std::endl;
    std::cout << "Plane(x=0) " << eval(pga3::e1) << " translated on the x-axis then meet with point on the x=0 plane"
              << eval(translated_plane ^ pga3::make_point(0, 3, 4)) << std::endl;

    auto t_z = pga3::translator(pga3::i, 2.0);
    auto r_45_x = pga3::rotor(pga3::e1, degrees2radians(45.)); // Rotate 45 around x axis
    auto r_90_x = pga3::rotor(pga3::e1, degrees2radians(90.)); // Rotate 90 around x axis
    auto rotated_plane = pga3::sandwich(pga3::e2, r_90_x);
    std::cout << "Plane (y=0) rotated around x-axis "<< eval(rotated_plane) << " then translated along z "
              << eval(pga3::sandwich(rotated_plane, t_z)) << std::endl;

}

void tests_with_lines()
{
    auto A = pga3::make_point(0, -1, 0);
    auto B = pga3::make_point(1, 1, -1);

    auto line_AB = pga3::line_from_points(A, B);
    print_line_info(line_AB, "line_AB");
    print_line_info(-1.0 * line_AB, "-line_AB");
    print_line_info(pga3::Line_t() - line_AB, "0-line_AB");

    auto line_object = pga3::Line(pga3::line_from_points(A, B));
    std::cout << "line_AB as object " << line_object << std::endl;

    std::cout << "e0 ^ line_AB " << eval(pga3::e0^line_AB) << std::endl;
    std::cout << "e1 ^ line_AB " << eval(pga3::e1^line_AB) << std::endl;
    std::cout << "e2 ^ line_AB " << eval(pga3::e2^line_AB) << std::endl;
    std::cout << "e3 ^ line_AB " << eval(pga3::e3^line_AB) << std::endl;

}

void test_bivectors()
{
    std::cout << "e12 " << eval(pga3::e12) << std::endl;
    std::cout << "e31 " << eval(pga3::e3 ^ pga3::e1) << std::endl;
    std::cout << "e23 " << eval(pga3::e23) << std::endl;

}


void test_trivectors()
{
    std::cout << "e123 " << eval(pga3::e123) << std::endl;
    std::cout << "e021 " << eval(pga3::e0 ^ pga3::e2 ^ pga3::e1) << std::endl;
    std::cout << "e013 " << eval(pga3::e013) << std::endl;
    std::cout << "e032 " << eval(pga3::e0 ^ pga3::e3 ^ pga3::e2) << std::endl;
}

pga3::space::algebra::element_t first_of_study_number(pga3::Line_t l)
{
    std::cout << "l: " << l  << std::endl;

    auto sn_B = pga3::bivector_axis(l);
    auto sn = sn_B.first.first *pga3::one + sn_B.first.second*pga3::I;
    auto sn_conj = sn_B.first.first * pga3::one - sn_B.first.second*pga3::I;
    std::cout << "Study number: " << sn_B.first.first << ", " << sn_B.first.second << " bv: " << eval(sn_B.second) << std::endl;
    std::cout << "Study number norm: " << ::sqrt(fabs((eval(sn * sn_conj)).template element<0>())) << std::endl;
    return sn_B.first.first;
}

template <typename TX, typename TY, typename TZ, typename DX, typename DY, typename DZ> inline
auto make_unnormalized_line(TX px, TY py, TZ pz, DX dx, DY dy, DZ dz)
{
    return (dx * pga3::dual_k + dy * pga3::dual_j + dz * pga3::dual_i + px * pga3::k + py * pga3::j + pz * pga3::i);
}

template<typename T1, typename T2, typename T3, typename T4>
inline
auto make_unnormalized_plane(T1 a, T2 b, T3 c, T4 d) {
    return d * pga3::e0 + a * pga3::e1 + b * pga3::e2 + c * pga3::e3;
}

void test_planes()
{

    auto pl = make_unnormalized_plane(7, 8, 9, 3);
    std::cout << "pl: " << pl << std::endl;
    std::cout << "pl: " << eval(pl) << std::endl;
    std::cout << "pl normalized: " << eval(pga3::normalize(pl)) << std::endl;
    std::cout << "Plane object pl: " << pga3::Plane(pl) << std::endl;
    std::cout << "Plane object pl normalized: " << pga3::Plane(pga3::Plane(pl).normalized()) << std::endl;

    auto ideal_pl = make_unnormalized_plane(0, 0, 0, 5);
    std::cout << "ideal_pl: " << ideal_pl << std::endl;
    std::cout << "ideal_pl: " << eval(ideal_pl) << std::endl;
    std::cout << "ideal_pl normalized: " << eval(pga3::normalize(ideal_pl)) << std::endl;
    std::cout << "Plane object ideal_pl: " << pga3::Plane(ideal_pl) << std::endl;
    std::cout << "Plane object ideal_pl normalized: " << pga3::Plane(pga3::Plane(ideal_pl).normalized()) << std::endl;
    std::cout << "ideal_pl normalized**2 " << eval(pga3::normalize(pga3::normalize(ideal_pl)) * pga3::normalize(pga3::normalize(ideal_pl))) << std::endl;

    auto ideal_pl2 = make_unnormalized_plane(0, 0, 0, -5);
    std::cout << "unnormalized ideal_pl2 " << eval(ideal_pl2) << std::endl;
    std::cout << "normalized ideal_pl2 " << eval(pga3::normalize(ideal_pl2)) << std::endl;
    std::cout << "normalized ideal_pl2**2 " << eval(pga3::normalize(ideal_pl2) * pga3::normalize(ideal_pl2)) << std::endl;

    auto pl2 = make_unnormalized_plane(3, 4, 5, 1);
    std::cout << "unnormalized pl2 " << eval(pl2) << std::endl;
    std::cout << "normalized pl2 " << eval(pga3::normalize(pl2)) << std::endl;
    std::cout << "normalized pl2**2 " << eval(pga3::normalize(pl2) * pga3::normalize(pl2)) << std::endl;

    auto pl3 = make_unnormalized_plane(3, 4, 5, 1);
    std::cout << "unnormalized pl3 " << eval(pl3) << std::endl;
    std::cout << "normalized pl3 " << eval(pga3::normalize(pl3)) << std::endl;
    std::cout << "normalized pl3**2 " << eval(pga3::normalize(pl3) * pga3::normalize(pl3)) << std::endl;
    std::cout << "unnormalized pl3**2 " << eval(pl3 * pl3) << std::endl;

    auto p = pga3::dual(pl3);
    std::cout << "unnormalized p " << eval(p) << std::endl;
    std::cout << "normalized p " << eval(pga3::normalize(p)) << std::endl;
    std::cout << "normalized p**2 " << eval(pga3::normalize(p) * pga3::normalize(p)) << std::endl;

}
void test_norms()
{
    auto P = pga3::make_point(2, 3, 4);
    auto Q = pga3::make_ideal_point(2, 3, 4);

    auto P_as_object = pga3::Point(P);
    std::cout << "P as object: " << P_as_object << " norm of P: " << pga3::norm(P) << std::endl;
    std::cout << "Q as object: " << pga3::Point(Q) << " norm of Q: " << pga3::norm(Q) << std::endl;
    std::cout << "normalized P: " << pga3::Point(pga3::Point(P).normalized()) << std::endl;
    std::cout << "normalized Q: " << pga3::Point(pga3::Point(Q).normalized()) << std::endl;
    std::cout << "normalized P as mv: " << eval(pga3::Point(P)) << std::endl;

    auto le = make_unnormalized_line(7, 8, 9, 0, 0, 0);
    auto li = make_unnormalized_line(0, 0, 0, 2, 3, 4);
    auto l = make_unnormalized_line(7, 8, 9, 2, 3, 4);

    std::cout << "line le: " << eval(le) << " norm of le: " << pga3::norm(le) << std::endl;
    std::cout << "line li: " << eval(li) << " norm of li: " << pga3::norm(li) << std::endl;
    std::cout << "line l: " << eval(l) << " norm of l: " << pga3::norm(l) << std::endl;

    std::cout << "Line object le: " << pga3::Line(le) << std::endl;
    std::cout << "Line object  li: " << pga3::Line(li) << std::endl;
    std::cout << "Line object  l: " << pga3::Line(l) << std::endl;
    std::cout << "Line object le normalized: " << pga3::Line(pga3::Line(le).normalized()) << std::endl;
    std::cout << "Line object  li: " << pga3::Line(pga3::Line(li).normalized()) << std::endl;
    std::cout << "Line object  l: " << pga3::Line(pga3::Line(l).normalized()) << std::endl;

    std::cout << "Point P + line l norm: " << pga3::norm(P + l) << std::endl;
    std::cout << "Point Q + line l norm: " << pga3::norm(Q + l) << std::endl;
    std::cout << "Point P + line le norm: " << pga3::norm(P + le) << std::endl;
    std::cout << "Point Q + line le norm: " << pga3::norm(Q + le) << std::endl;
    std::cout << "Point P + line li norm: " << pga3::norm(P + li) << std::endl;
    std::cout << "Point Q + line li norm: " << pga3::norm(Q + li) << std::endl;
    std::cout << "line le ^ line li: " << (le ^ li) << std::endl;
    std::cout << "line le ^ line li norm: " << pga3::norm(le ^ li) << std::endl;
    std::cout << "line le V line li: " << pga3::vee(le, li) << std::endl;
    std::cout << "line le V line li norm: " << pga3::norm(pga3::vee(le, li)) << std::endl;

    auto a = pga3::normalize(make_unnormalized_plane(0, 3, 0, 0));
    std::cout << "plane a: " << a << std::endl;
    std::cout << "a norm: " << pga3::norm(a) << std::endl;
    std::cout << "a*a: " << (a*a) << std::endl;

    auto A = pga3::make_point(0, 0, 3);
    std::cout << "Point A: " << A << std::endl;
    std::cout << "Polar Point A: " << (A*pga3::I) << std::endl;
    std::cout << "A norm: " << pga3::norm(A) << std::endl;
    std::cout << "A normalized: " << pga3::normalize(A) << std::endl;
    std::cout << "pga3::normalize(A)**2: " << (pga3::normalize(A)*pga3::normalize(A)) << std::endl;
    auto B = pga3::make_point(3, 0, 3);
    std::cout << "Point B: " << B << std::endl;
    std::cout << "Polar Point B: " << (B*pga3::I) << std::endl;
    auto C = pga3::make_point(3, 4, 3);
    auto AB = pga3::line_from_points(A, B);
    auto BC = pga3::line_from_points(B, C);
    auto CA = pga3::line_from_points(C, A);

    std::cout << "Length of line AB: " << pga3::norm(AB) << std::endl;
    std::cout << "Length of line BC: " << pga3::norm(BC) << std::endl;
    std::cout << "Length of line CA: " << pga3::norm(CA) << std::endl;
    std::cout << "Area of triangle ABC: " << 0.5 * pga3::norm(AB + BC + CA) << std::endl;

    auto surface_element = (a + A);
    auto polar_surface_element = surface_element * pga3::I;
    std::cout << "a ^ A: " << (a ^ A) << std::endl;
    std::cout << "Surface element a + A: " << surface_element << std::endl;
    std::cout << "Polar surface element aI + AI: " << polar_surface_element << std::endl;
    std::cout << "Norm of surface element a + A: " << pga3::norm(surface_element) << std::endl;
    std::cout << "Norm of polar surface element aI + AI: " << pga3::norm(polar_surface_element) << std::endl;
    std::cout << "(a + A)**2: " << (surface_element*surface_element) << std::endl;
    std::cout << "(aI + AI)**2: " << (polar_surface_element*polar_surface_element) << std::endl;
    auto axis = (a ^ A*pga3::I);
    std::cout << "axis: a ^ A*I: " << axis << std::endl;
    auto spear = pga3::vee(A, a*pga3::I);
    std::cout << "spear: A V a*I: " << spear << std::endl;
    std::cout << "a ^ axis: " << (a ^ axis) << std::endl;
    std::cout << "A V spear: " << pga3::vee(A, spear) << std::endl;

    std::cout << "dual of (1e0+x*1e1+y*1e2+z*1e3)(3,4,5)" <<
                 pga3::dual((pga3::e0+3*pga3::e1+4*pga3::e2 + 5*pga3::e3)) << " make_point: " <<
                 pga3::make_point(3,4,5) << std::endl;

//    std::cout << "Norm of -3 : " << pga3::norm(-3.0 * pga3::one ) << std::endl;
//    std::cout << "Scalar - 3 * pseudoscalar: " << pga3::norm(pga3::one - 3.0 * pga3::I) << std::endl;
//    std::cout << "3 * pseudoscalar: " << pga3::norm(3.0 * pga3::I) << std::endl;


}

int main()
{
//    print_cayley_table();
//    print_info_all_basis();
//    print_bivector_info();
//    print_plane_info();
//    test_with_points();
//    tests_with_lines();
//    test_planes();
//    test_trivectors();
    test_norms();

}
