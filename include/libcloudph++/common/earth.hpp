#pragma once

#include <libcloudph++/common/units.hpp>
#include <libcloudph++/common/macros.hpp>
#include <libcloudph++/common/moist_air.hpp>//TODO better name

namespace libcloudphxx
{
  namespace common
  {
    namespace earth
    {
      using moist_air::R_d;

      // acceleration due to gravity
      libcloudphxx_const(si::acceleration, g, 9.81, si::metres_per_second_squared)

      //standard pressure and temperature (ICAO)
      libcloudphxx_const(si::pressure, p_stp, 101325, si::pascal)
      libcloudphxx_const(si::temperature, T_stp, 273.15 + 15, si::kelvins)

      libcloudphxx_const_derived(si::mass_density, rho_stp, p_stp<real_t>() / T_stp<real_t>() / R_d<real_t>())
    };
  };
};
