#include "gaalet.h"
#include "pga3.h"
#include <cmath>
#include <vector>

using namespace gaalet;


template <class MV>
void print_basis_info(const MV& p, const std::string& name) {
    std::cout << name << ": " << p 
//              << " polar of " << name << ": " << pga3::polar(p) 
              << " dual of " << name << ": " << ::dual(p) 
              << " square of " << name << ": " << p * p // + 0. * pga3::one)
              << std::endl;
}


void print_point_info(const pga3::Point_t& p, const std::string& name) {
    std::cout << name << ": " << p << " polar of " << name << ": " << pga3::polar(p) << " dual of " << name << ": " << ::dual(p) << std::endl;
    std::cout << "normalized " << name << ": " << normalize(p) << std::endl;
    std::cout << name << "[0]: " << p[0] << " " << name << "[1]: " << p[1] << " " 
              << name << "[2]: " << p[2] << " " << name << "[3]: " << p[3] << std::endl;
}


void print_line_info(const pga3::Line_t& l, const std::string& name) {
    std::cout << name << ": " << l << " polar of " << name << ": " << pga3::polar(l) << " dual of " << name << ": " << ::dual(l) << std::endl;
    std::cout << "normalized " << name << ": " << normalize(l) << std::endl;
    std::cout << name << "[0]: " << l[0] << " " << name << "[1]: " << l[1] << " " << name << "[2]: " << l[2] << std::endl;
    std::cout << name << "[3]: " << l[3] << " " << name << "[4]: " << l[4] << " " << name << "[5]: " << l[5] << std::endl;
}


int main()
{
//    std::cout << pga3::e0 << std::endl;

    std::cout << pga3::space::algebra::metric::dimension << std::endl;
    std::cout << pga3::space::algebra::metric::signature_bitmap << std::endl;
    std::cout << pga3::space::algebra::metric::degenerate_bitmap << std::endl;
//
    print_basis_info(pga3::one, "1");
    print_basis_info(pga3::e0, "e0");
    print_basis_info(pga3::e1, "e1");
    print_basis_info(pga3::e2, "e2");
    print_basis_info(pga3::e3, "e3");
    print_basis_info(pga3::e01, "e01");
    print_basis_info(pga3::e02, "e02");
    print_basis_info(pga3::e03, "e03");
    print_basis_info(pga3::e12, "e12");
    print_basis_info(pga3::e13, "e13");
    print_basis_info(pga3::e23, "e23");
    print_basis_info(pga3::e012, "e012");
    print_basis_info(pga3::e013, "e013");
    print_basis_info(pga3::e023, "e023");
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
    print_basis_info(pga3::EX, "EX");
    print_basis_info(pga3::EY, "EY");
    print_basis_info(pga3::EZ, "EZ");
    
    std::cout << "::polar(pga3::one) " << ::polar(pga3::one) << std::endl;
    std::cout << "::polar(pga3::I) " << ::polar(pga3::I) << std::endl;
    std::cout << "::dual(pga3::one) " << ::dual(pga3::one) << std::endl;
    std::cout << "::dual(pga3::I) " << ::dual(pga3::I) << "::dual(::dual(pga3::I)) " << ::dual(::dual(pga3::I)) << std::endl;
    std::cout << "::dual(pga3::I) " << ::dual(pga3::I) << std::endl;
    std::cout << ~pga3::one << std::endl;
    std::cout << ~pga3::I << std::endl;
    std::cout << " " << std::endl;
    
    std::cout << (pga3::e1 ^ pga3::e2 ^ pga3::e0) << std::endl;

//    auto a_line = pga3::line(2, 3, 4, 5, 6, 7);
//    auto a_plane = pga3::plane(2, 3, 4, 7);
//    auto a_point = pga3::point(2, 3, 4);
//    auto an_ideal_point = pga3::ideal_point(2, 3, 4);
//    auto another_point = pga3::point(4, 7, 8) * 2.;
//    auto A = pga3::point(0, -1, 0), B = pga3::point(1, 1, -1), C = pga3::point(-1, 1, -1), D = pga3::point(1, 1, 1),
//         E = pga3::point(-1, 1, 1), F = pga3::point(5, 6, 7);
//    auto centroid = A + B + C + D + E;
//
//    std::cout << "polar(e1 ^ e2): " << pga3::polar(pga3::e1 ^ pga3::e2) << std::endl;
//    std::cout << "dual(e1 ^ e2): " << ::dual(pga3::e1 ^ pga3::e2) << std::endl;
//    std::cout << "polar(e1): " << pga3::polar(pga3::e1) << std::endl;
//    std::cout << "polar(E1): " << pga3::polar(pga3::E1) << std::endl;
//    std::cout << "dual(e1): " << ::dual(pga3::e1) << std::endl;
//    std::cout << "dual(E1): " << ::dual(pga3::E1) << std::endl;
//    std::cout << "dual(1): " << ::dual(pga3::one) << std::endl;
//    std::cout << "dual(I): " << ::dual(pga3::I) << std::endl;
//    std::cout << "polar(1): " << pga3::polar(pga3::one) << std::endl;
//    std::cout << "polar(I): " << pga3::polar(pga3::I) << std::endl;
//
//    std::cout << "I " << (pga3::I) << std::endl;
//    std::cout << "polar(e012)" << pga3::polar(pga3::e012) << std::endl;
//    std::cout << "e013 * I " << pga3::e013 * pga3::I << std::endl;
//    std::cout << "e023 * I " << pga3::e023 * pga3::I << std::endl;
//    std::cout << "e123 * I " << pga3::e123 * pga3::I << std::endl;
//    std::cout << "polar(e123) " << pga3::polar(pga3::e123) << std::endl;
//
//    print_point_info(F, "F");
//    print_point_info(A, "A");
//    print_point_info(B, "B");
//    std::cout << "(A ^ B): " << (A ^ B) << std::endl;
//    std::cout << "(polar(A) ^ polar(B)): " << (pga3::polar(A) ^ pga3::polar(B)) << std::endl;
//    std::cout << "vee(A, B): " << pga3::vee(A, B) << std::endl;
//
//    print_point_info(centroid, "centroid");
//    std::cout << "inverse centroid: " << ::inverse(normalize(centroid)) << std::endl;
//    
//    std::cout << "A line from points(A, B): " << pga3::line_from_points(A, B) << std::endl;
//    print_line_info(pga3::line_from_points(A, B), "line");
//
//
//    std::cout << "A line from points(centroid,C+E): " << pga3::line_from_points(centroid, C + E) << std::endl;
//    std::cout << "Sum of points: " << normalize(a_point + another_point) << std::endl;
//    std::cout << "A point doubled then normalized: " << normalize(2. * a_point) << std::endl;
//    std::cout << "A point normalized: " << normalize(a_point) << std::endl;
//    std::cout << "Another point normalized: " << normalize(another_point) << std::endl;
//    std::cout << "An ideal point dualized: " << ::dual(an_ideal_point) << std::endl;
//    // Next line fails to compile since can't normalize an ideal point
//    //    std::cout << "An ideal point normalized: " << normalize(an_ideal_point) << std::endl;
//    std::cout << "A point dualized: " << ::dual(a_point) << std::endl;
//    std::cout << "A line normalized: " << normalize(a_line) << std::endl;
//    std::cout << "A line dualized: " << ::dual(a_line) << std::endl;
//    std::cout << "A plane normalized: " << normalize(a_plane) << std::endl;
//    std::cout << "A plane dualized, then normalized: " << normalize(::dual(a_plane)) << std::endl;
//    std::cout << "A line from points: " << pga3::line_from_points(a_point, another_point) << std::endl;
//    std::cout << "Squared line from points: "
//              << pga3::line_from_points(a_point, another_point) * pga3::line_from_points(a_point, another_point)
//              << std::endl;
//    std::cout << "Squared line from points: "
//              << pga3::line_from_points(a_point, another_point) * pga3::line_from_points(a_point, another_point)
//              << std::endl;
//    std::cout << "A plane from points: " << pga3::plane_from_points(A, B, C) << std::endl;
//    std::cout << "A line from planes: "
//              << pga3::line_from_planes(pga3::plane_from_points(A, B, C), pga3::plane_from_points(B, C, D))
//              << std::endl;
//    std::cout << "A point from line and plane: "
//              << pga3::point_from_line_and_plane(
//                     pga3::line_from_points(centroid, C + E), pga3::plane_from_points(A, D, B))
//              << std::endl;
//    std::cout << "sum of lines:" << (pga3::line_from_points(B, D) + pga3::line_from_points(C, E)) << std::endl;
//    std::cout << "motor sandwich applied to a line:"
//              << pga3::sandwich(
//                     pga3::line_from_points(C, E), pga3::motor(pga3::line_from_points(centroid, C + E), 2, M_PI / 3.))
//              << std::endl;
////    std::cout << "motor2 sandwich applied to a line:"
////              << pga3::sandwich(
////                     pga3::line_from_points(C, E), pga3::screw2(pga3::line_from_points(centroid, C + E), 2, M_PI / 3.))
////              << std::endl;
//    std::cout << "motor sandwich applied to point:" << A << " "
//              << (pga3::sandwich(A, pga3::motor(pga3::line_from_points(centroid, C + E), 2, M_PI / 3.))) << std::endl;
////    std::cout << "motor2 sandwich applied to point:" << A << " "
////              << (pga3::sandwich(A, pga3::screw2(pga3::line_from_points(centroid, C + E), 2, M_PI / 3.))) << std::endl;
//    std::cout << "translator sandwich applied to point:" << A << " "
//              << (pga3::sandwich(A, pga3::translator(pga3::line_from_points(centroid, C + E), 1))) << std::endl;
//    std::cout << "vee applied to line:" << pga3::vee(pga3::E0, 0.5 * normalize(pga3::line_from_points(centroid, C + E)))
//              << std::endl;
//
//        std::cout << "Pseudoscalar normalized: " << normalize(-3*pga3::I) << std::endl;
}
