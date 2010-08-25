#include "gaalet.h"

int main()
{
   typedef typename gaalet::insert_element_to_melist<1, 2, gaalet::mel_null, -1>::melist melist;

   std::cout << "size clist: " << melist::clist::size << std::endl;
   std::cout << "element 0: " << get_element<0, typename melist::clist>::value << std::endl;
}
