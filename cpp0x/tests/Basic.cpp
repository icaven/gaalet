#include "gaalet.h"

template <class MV>
void print_basis_info(const MV& p, const std::string& name) {
    std::cout << name << ": " << p 
              << " dual of " << name << ": " << ::dual(p) 
              << " square of " << name << ": " << eval(p * p )
              << std::endl;
}

// Convenience function and macro to generate the configuration list values
inline constexpr gaalet::conf_t conf_bits_e(const int v=-1) { return v < 0 ? 0 : 1 << (v-1); };
#define conf_bits(v1, v2) (conf_bits_e(v1)  | conf_bits_e(v2) )

int main()
{
   gaalet::mv<3, 5>::type a = {1.0, 2.0};
   gaalet::mv<3, 4>::type b = {4.0, 7.0};

   auto c = a - b;

   std::cout << "a: " << a << std::endl;
   std::cout << "b: " << b << std::endl;
   std::cout << "c: " << c << std::endl;

   gaalet::mv<0>::type d;

   typedef gaalet::algebra<gaalet::signature<3,0> > em;
   em::mv<1,2,4>::type e = {1,2,3};
   em::mv<1,2,4>::type f = {3,4,5};
   std::cout << "ef: " << e*f << std::endl;
   
   
    // The basis elements of the algebra

    // Scalar
    const em::mv<0x0>::type one = { 1.0 };

    // Vectors
    const em::mv<conf_bits_e(1)>::type e1 = { 1.0 };
    const em::mv<conf_bits_e(2)>::type e2 = { 1.0 };
    const em::mv<conf_bits_e(3)>::type e3 = { 1.0 };

    // Bivectors
    const em::mv<conf_bits(1, 2)>::type e12 = (e1 ^ e2);
    const em::mv<conf_bits(1, 3)>::type e13 = (e1 ^ e3);
    const em::mv<conf_bits(2, 3)>::type e23 = (e2 ^ e3);

    // The pseudoscalar (representing all space)
    const gaalet::conf_t pseudoscalar_conf = gaalet::Power<2, em::algebra::metric::dimension>::value-1;
    const em::mv<pseudoscalar_conf>::type I = (e1 ^ e2 ^ e3);
   
    print_basis_info(one, "1");
    print_basis_info(e1, "e1");
    print_basis_info(e2, "e2");
    print_basis_info(e3, "e3");
    print_basis_info(e12, "e12");
    print_basis_info(e13, "e13");
    print_basis_info(e23, "e23");
    print_basis_info(I, "I");

}
