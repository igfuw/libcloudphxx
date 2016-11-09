//#include <libcloudph++/common/detail/fake_units.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/pow.hpp>
#include <boost/units/systems/si.hpp>

#include <cassert>

#include <iostream>

//namespace si = libcloudphxx::common::detail::fake_units::si;
//using libcloudphxx::common::detail::fake_units::quantity;
//using libcloudphxx::common::detail::fake_units::pow;

using namespace boost::units;

template <typename real_t>
void test()
{
  // empty constructor
  quantity<si::pressure, real_t> p0;
  quantity<power_typeof_helper<si::pressure, static_rational<2>>::type, real_t> p00;

  // ctor with assignment 
  quantity<si::pressure, real_t> p1 = static_cast<real_t>(200) * si::pascals;
  if (p1.value() != static_cast<real_t>(200)) assert(false),throw;

  // assignment
  p0 = p1;
  if (p0.value() != p1.value()) assert(false),throw;

  // multiplication
  p00 = p1 * p1;
  if (p00.value() != p1.value() * p1.value()) assert(false),throw;

  // division
  p0 = p00 / p1;
  if (p0.value() != p1.value()) assert(false),throw;
  
  // power
  p0 = boost::units::pow<1>(p1);
  if (p0.value() != 200) assert(false),throw;

  // addition
  p0 = p1 + p1;
  if (p0.value() != 400) assert(false),throw;

  // subtraction
  p0 = p1 - p1;
  if (p0.value() != 0) assert(false),throw;

  // cast
  double a = static_cast<double>(p0.value());
  if (a != p0.value()) assert(false),throw;

  // int over quantity
  p0 = p00 * (static_cast<real_t>(1) / p1);
  if (p0.value() != 200) assert(false),throw;

  // int times quantity
  p00 = p1 * (static_cast<real_t>(2) * p1);
  if (p00.value() != 200*(2. * 200)) assert(false),throw;

  // q / int
  p00 = p1 * (p1 / static_cast<real_t>(2));
  if (p00.value() != 200*(200. / 2)) assert(false),throw;

  // q * int
  p00 = p1 * (p1 * static_cast<real_t>(2));
  if (p00.value() != 200*(2. * 200)) assert(false),throw;

/*
  // addition of int
  p0 = p1 * (p1 + 1);
  if (p0.value() != 200*201) assert(false),throw;

  // addition to int
  p0 = p1 * (static_cast<real_t>(1) + p1);
  if (p0.value() != 200*201) assert(false),throw;

  // int minus quantity
  p0 = p1 * (static_cast<real_t>(2) - p1);
  if (p0.value() != 200*(2. - 200)) assert(false),throw;

  // q + int
  p0 = p1 * (p1 + static_cast<real_t>(1));
  if (p0.value() != 200*201) assert(false),throw;

  // quantity - int
  p0 = p1 * (p1 - static_cast<real_t>(2));
  if (p0.value() != 200*(200-2)) assert(false),throw;
*/
}

int main()
{ 
  test<float>();
  test<double>();
  test<long double>();
}
