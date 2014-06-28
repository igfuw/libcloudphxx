#pragma once

#include <libcloudph++/common/units.hpp>
#include <libcloudph++/common/macros.hpp> 

namespace libcloudphxx
{
  namespace common
  {
    typedef divide_typeof_helper<si::mass, si::amount>::type mass_over_amount;

    namespace molar_mass
    {
      libcloudphxx_const(mass_over_amount, M_H,     1*1e-3, si::kilograms / si::moles)
      libcloudphxx_const(mass_over_amount, M_OH,   17*1e-3, si::kilograms / si::moles)
      libcloudphxx_const(mass_over_amount, M_SO2,  64*1e-3, si::kilograms / si::moles)
      libcloudphxx_const(mass_over_amount, M_H2O2, 34*1e-3, si::kilograms / si::moles)
      libcloudphxx_const(mass_over_amount, M_O3,   48*1e-3, si::kilograms / si::moles)
    };
  };
};
