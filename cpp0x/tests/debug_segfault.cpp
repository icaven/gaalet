#include "gaalet.h"

// This program was used to debug the reference overwrite problem
// and is now just for documentation purposes

typedef gaalet::algebra<gaalet::signature<3,0,1>> ma;

int main()
{
   ma::mv<0>::type one={1.0};
   ma::mv<1>::type e1={1.0};
   ma::mv<2>::type e2={1.0};
   ma::mv<4>::type e3={1.0};
   ma::mv<8>::type e0={1.0};
   ma::mv<0xf>::type I= (e1^e2^e3^e0);


   auto x = 11.0*e2*e3 + 12.0*e3*e1 + 13.0*e1*e2;
   // x will be an expression object structured in the form of a tree with left and right branches
   //             add
   //          /      \
   //       add         gp
   //       /  \       /  \
   //    gp     gp    svp  mv
   //   /  \    / \   / \  |
   // svp  mv  svp mv 13 e1 e2
   // / \  |    /\  |
   //11 e2 e3 12 e3 e1
   // where: add is addition, gp is geometric product, svp is scalar multivector product, mv is multivector
   //
   
   // After the I*x operation on the next line, parts of the x object have been overwritten
   // but they were referencing other parts of x
   auto X = one + I*x;
   std::cout << "X: " << X << std::endl;
   std::cout << std::endl;
}
