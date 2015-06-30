#include <libcloudph++/common/detail/toms748.hpp>

#include <iostream>

double f(double x)
{
  return x*x - 2;
}

int main()
{
  double a=1, b=2;
  uintmax_t maxiter=10;
  auto pair = libcloudphxx::common::detail::toms748_solve(
      f, 
      a, 
      b, 
      libcloudphxx::common::detail::eps_tolerance<double>(sizeof(double) * 8 / 2), 
      maxiter
    ); 
  std::cerr << (pair.first + pair.second)/2 << std::endl;
}
