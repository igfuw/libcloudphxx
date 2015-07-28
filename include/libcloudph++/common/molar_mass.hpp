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
      //trace gases
      libcloudphxx_const(mass_over_amount, M_SO2,  64*1e-3, si::kilograms / si::moles)
      libcloudphxx_const(mass_over_amount, M_H2O2, 34*1e-3, si::kilograms / si::moles)
      libcloudphxx_const(mass_over_amount, M_O3,   48*1e-3, si::kilograms / si::moles)

      // dissolved trace gases
      libcloudphxx_const(mass_over_amount, M_H2O2_H2O, (34+18) * 1e-3, si::kilograms / si::moles)
      libcloudphxx_const(mass_over_amount, M_O3_H2O,   (48+18) * 1e-3, si::kilograms / si::moles)

      //H2O
      libcloudphxx_const(mass_over_amount, M_H,     1*1e-3, si::kilograms / si::moles)
      libcloudphxx_const(mass_over_amount, M_OH,   17*1e-3, si::kilograms / si::moles)
  
      //SO2 * H2O
      libcloudphxx_const(mass_over_amount, M_SO2_H2O, 82*1e-3, si::kilograms / si::moles)
      libcloudphxx_const(mass_over_amount, M_HSO3,    81*1e-3, si::kilograms / si::moles)
      libcloudphxx_const(mass_over_amount, M_SO3,     80*1e-3, si::kilograms / si::moles)

      //H2SO4
      libcloudphxx_const(mass_over_amount, M_H2SO4, 98*1e-3, si::kilograms / si::moles) 
      libcloudphxx_const(mass_over_amount, M_HSO4,  97*1e-3, si::kilograms / si::moles)
      libcloudphxx_const(mass_over_amount, M_SO4,   96*1e-3, si::kilograms / si::moles)
    };
  };
};
