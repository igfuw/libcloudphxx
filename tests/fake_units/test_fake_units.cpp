#include <libcloudph++/common/detail/fake_units.hpp>
#include <cassert>

namespace si = libcloudphxx::common::detail::fake_units::si;
using libcloudphxx::common::detail::fake_units::quantity;
using libcloudphxx::common::detail::fake_units::pow;

int main()
{
  // empty constructor
  quantity<si::pressure, double> p0;

  // ctor with assignment 
  quantity<si::pressure, double> p1 = double(200) * si::pascals;
  assert(p1.value == double(200));

  // assignment
  p0 = p1;
  assert(p0.value == p1.value);

  // multiplication
  p0 = p1 * p1;
  assert(p0.value == p1.value * p1.value);

  // division
  p0 = p1 / p1;
  assert(p0.value == 1);
  
  // power
  p0 = pow(p1, p0);
  assert(p0.value == 200);

  // addition
  p0 = p1 + p1;
  assert(p0.value == 400);

  // conversion
  double a = p0;
}
