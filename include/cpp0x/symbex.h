#ifndef __GAALET_SYMBOL_H
#define __GAALET_SYMBOL_H

#include "string"
#include "sstream"

namespace gaalet
{

///symbex: symbolic expression
struct symbex {
   symbex(const std::string& expr_ = "null")
      :  expr(expr_)
   { }

   //template<typename A>
   //symbex(const A& expr_)
   symbex(const double& expr_)
   {
      std::stringstream as;
      as << expr_;
      expr = as.str();
   }

   //template<typename A>
   //symbex(const A& expr_)
   symbex(const char expr_[])
   {
      std::stringstream as;
      as << expr_;
      expr = as.str();
   }

   std::string expr;
};

} //end namespace gaalet

std::ostream& operator<<(std::ostream& os, const gaalet::symbex& s)
{
   os << s.expr;
   return os;
}

gaalet::symbex
operator+(const gaalet::symbex& l, const gaalet::symbex& r) {
   return gaalet::symbex("(" + l.expr + "+" + r.expr + ")");
}

gaalet::symbex
operator-(const gaalet::symbex& l, const gaalet::symbex& r) {
   return gaalet::symbex("(" + l.expr + "-" + r.expr + ")");
}


gaalet::symbex
operator*(const gaalet::symbex& l, const gaalet::symbex& r) {
   return gaalet::symbex(l.expr + "*" + r.expr);
}

gaalet::symbex
operator*(const double& l, const gaalet::symbex& r) {
   if(l == 1.0) {
      return gaalet::symbex(r.expr);
   }
   else if(l == -1.0) {
      return gaalet::symbex("(-" + r.expr+")");
   }
   else if(l == 0.0) {
      return gaalet::symbex("");
   }
   else {
      std::stringstream ls;
      ls << l;

      return gaalet::symbex(ls.str() + "*" + r.expr);
   }
}

gaalet::symbex
operator*(const gaalet::symbex& l, const double& r) {
   if(r == 1.0) {
      return gaalet::symbex(l.expr);
   }
   else if(r == -1.0) {
      return gaalet::symbex("(-" + l.expr+")");
   }
   else if(r == 0.0) {
      return gaalet::symbex("");
   }
   else {
      std::stringstream rs;
      rs << r;

      return gaalet::symbex(l.expr + "*" + rs.str());
   }
}


gaalet::symbex
operator/(const gaalet::symbex& l, const gaalet::symbex& r) {
   return gaalet::symbex(l.expr + "/" + r.expr);
}


gaalet::symbex
sin(const gaalet::symbex& a) {
   return gaalet::symbex("sin(" + a.expr + ")");
}

gaalet::symbex
cos(const gaalet::symbex& a) {
   return gaalet::symbex("cos(" + a.expr + ")");
}

gaalet::symbex
sqrt(const gaalet::symbex& a) {
   return gaalet::symbex("sqrt(" + a.expr + ")");
}

#endif
