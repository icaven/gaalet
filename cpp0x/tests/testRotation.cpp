#include "gaalet.h"

int main()
{
   auto a = 1*e1 + 2*e2 + 3*e3;
   auto n = 0.1*e1 + 0.2*e2 + 0.3*e3;
   auto R = exp( (-0.5)*i*n );
   auto b = eval(grade<1>(R*a*(!R)));
}
