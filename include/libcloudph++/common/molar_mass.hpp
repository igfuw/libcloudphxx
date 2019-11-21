#pragma once

#include "units.hpp"
#include "macros.hpp" 

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
      libcloudphxx_const(mass_over_amount, M_NH3,  17*1e-3, si::kilograms / si::moles)
      libcloudphxx_const(mass_over_amount, M_HNO3, 63*1e-3, si::kilograms / si::moles)
      libcloudphxx_const(mass_over_amount, M_CO2,  44*1e-3, si::kilograms / si::moles)
                                                  //sic!
      //H2O
      libcloudphxx_const(mass_over_amount, M_H,     1*1e-3, si::kilograms / si::moles)
      libcloudphxx_const(mass_over_amount, M_OH,   17*1e-3, si::kilograms / si::moles)
  
      //SO2 * H2O
      libcloudphxx_const(mass_over_amount, M_SO2_H2O, 82*1e-3, si::kilograms / si::moles)
      libcloudphxx_const(mass_over_amount, M_HSO3,    81*1e-3, si::kilograms / si::moles)
      libcloudphxx_const(mass_over_amount, M_SO3,     80*1e-3, si::kilograms / si::moles)

      //NH3 * H2O
      libcloudphxx_const(mass_over_amount, M_NH3_H2O, 35*1e-3, si::kilograms / si::moles)
      libcloudphxx_const(mass_over_amount, M_NH4,     18*1e-3, si::kilograms / si::moles)

      //HNO3(aq)
      libcloudphxx_const(mass_over_amount, M_NO3,     62*1e-3, si::kilograms / si::moles)

      //CO2 * H2O
      libcloudphxx_const(mass_over_amount, M_CO2_H2O, 62*1e-3, si::kilograms / si::moles)
      libcloudphxx_const(mass_over_amount, M_HCO3,    61*1e-3, si::kilograms / si::moles)
      libcloudphxx_const(mass_over_amount, M_CO3,     60*1e-3, si::kilograms / si::moles)

      //H2SO4
      libcloudphxx_const(mass_over_amount, M_H2SO4, 98*1e-3, si::kilograms / si::moles) 
      libcloudphxx_const(mass_over_amount, M_HSO4,  97*1e-3, si::kilograms / si::moles)
      libcloudphxx_const(mass_over_amount, M_SO4,   96*1e-3, si::kilograms / si::moles)
    };
  };
};
