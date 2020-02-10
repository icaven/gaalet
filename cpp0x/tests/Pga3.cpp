
#include <cmath>
#include <vector>

#include "pga3_streaming.h"

#include "gaalet.h"
#include "pga3.h"
#include "pga3_ops.h"
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

void test_norms()
{
    auto P = pga3::make_point(2, 3, 4);
    auto Q = pga3::make_ideal_point(2, 3, 4);

    auto P_as_object = pga3::Point(P);
    std::cout << "P as object: " << P_as_object << std::endl;
    std::cout << "Q as object: " << pga3::Point(Q) << std::endl;
    std::cout << "normalized P: " << pga3::Point(pga3::Point(P).normalized()) << std::endl;
    std::cout << "normalized Q: " << pga3::Point(pga3::Point(Q).normalized()) << std::endl;
    std::cout << "normalized P as mv: " << eval(pga3::Point(P)) << std::endl;

    auto le = make_unnormalized_line(7, 8, 9, 0, 0, 0);
    auto li = make_unnormalized_line(0, 0, 0, 2, 3, 4);
    auto l = make_unnormalized_line(7, 8, 9, 2, 3, 4);

    std::cout << "line le: " << eval(le) << std::endl;
    std::cout << "line li: " << eval(li) << std::endl;
    std::cout << "line l: " << eval(l) << std::endl;

    std::cout << "Line object le: " << pga3::Line(le) << std::endl;
    std::cout << "Line object  li: " << pga3::Line(li) << std::endl;
    std::cout << "Line object  l: " << pga3::Line(l) << std::endl;
    std::cout << "Line object le normalized: " << pga3::Line(pga3::Line(le).normalized()) << std::endl;
    std::cout << "Line object  li: " << pga3::Line(pga3::Line(li).normalized()) << std::endl;
    std::cout << "Line object  l: " << pga3::Line(pga3::Line(l).normalized()) << std::endl;

    auto pl = make_unnormalized_plane(7, 8, 9, 1);
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

}

int main()
{
    print_cayley_table();
    print_info_all_basis();
    print_bivector_info();
    print_plane_info();
    test_with_points();
    tests_with_lines();
    test_trivectors();
    test_norms();

}
