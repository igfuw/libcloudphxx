#include <libcloudph++/common/detail/fake_units.hpp>
#include <cassert>

namespace si = libcloudphxx::common::detail::fake_units::si;
using libcloudphxx::common::detail::fake_units::quantity;
using libcloudphxx::common::detail::fake_units::pow;

template <typename real_t>
void test()
{
  // empty constructor
  quantity<si::pressure, double> p0;

  // ctor with assignment 
  quantity<si::pressure, double> p1 = double(200) * si::pascals;
  if (p1.value != double(200)) throw;

  // assignment
  p0 = p1;
  if (p0.value != p1.value) throw;

  // multiplication
  p0 = p1 * p1;
  if (p0.value != p1.value * p1.value) throw;

  // division
  p0 = p1 / p1;
  if (p0.value != 1) throw;
  
  // power
  p0 = pow(p1, p0);
  if (p0.value != 200) throw;

  // addition
  p0 = p1 + p1;
  if (p0.value != 400) throw;

  // subtraction
  p0 = p1 - p1;
  if (p0.value != 0) throw;

  // cast
  double a = p0;
  if (a != p0.value) throw;
}

int main()
{ 
  test<float>();
  test<double>();
  test<long double>();
}
