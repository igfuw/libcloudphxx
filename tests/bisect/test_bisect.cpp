#include <libcloudph++/common/detail/bisect.hpp>

#include <iostream>

double f(double x)
{
  return x*x - 2;
}

int main()
{
  double a=1, b=2, tol=.00001;
  std::cerr 
    << libcloudphxx::common::detail::bisect(f, a, b, tol) 
//    << " "
//    << libcloudphxx::common::detail::secant(f, a, b, tol) 
    << std::endl;
}
