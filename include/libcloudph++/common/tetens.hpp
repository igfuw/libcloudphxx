// Approximate saturation vapour pressure and saturation mixing ratio using the Tetens formula

#pragma once

#include <libcloudph++/common/units.hpp> // TODO: do detail?

namespace libcloudphxx
{
  namespace common
  {
    namespace tetens
    {
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::pressure, real_t> p_vs(
        const quantity<si::temperature, real_t> &T
      ) 
      {
        const real_t T_in_celsius(T / si::kelvins - 273.15);
        assert(T_in_celsius > 0.);
	return 6.1078e3 * exp( (17.27 * T_in_celsius) / ( T_in_celsius + 237.3)) * si::pascals;
      }

      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::dimensionless, real_t> r_vs(
        const quantity<si::temperature, real_t> &T,
        const quantity<si::pressure, real_t> &p // total pressure
      ) 
      {
        const real_t T_in_celsius(T / si::kelvins - 273.15);
        assert(T_in_celsius > 0.);
        // r_vs=3.8/(p*exp(-17.2693882*(T-273.15)/(T-35.86))-6.109)  p in mb, T in Kelvins
	return 380. / (p / si::pascals * exp(-17.2693882 * (T_in_celsius) / (T / si::kelvins - 35.86)) - 610.9);
      }
    };
  };
};
