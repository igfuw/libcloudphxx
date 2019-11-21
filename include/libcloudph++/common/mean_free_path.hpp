#pragma once

#include "units.hpp"
#include "macros.hpp"
#include "moist_air.hpp"

namespace libcloudphxx
{
  namespace common
  {
    namespace mean_free_path
    {
      // eq. 6.6 in Williams and Loyalka 1991
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::length, real_t> lambda_D(
        quantity<si::temperature, real_t> T  // temperature
      ) 
      {
        using moist_air::R_v;
        using moist_air::D_0; // TODO: perhaps D as a parameter to allow it to vary with T,p?
using std::pow;

        return real_t(2) * D_0<real_t>() / (
          pow(
            real_t(2) * R_v<real_t>() * T
            / si::metres_per_second / si::metres_per_second, // boost::pow works for double only :(
            real_t(.5)
          ) * si::metres_per_second
        );
      }

      // eq. 6.33 in Williams and Loyalka 1991
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::length, real_t> lambda_K(
        const quantity<si::temperature, real_t> &T,  // temperature
        const quantity<si::pressure, real_t> &p   // pressure
      ) 
      {
        using moist_air::R_d;
        using moist_air::K_0;
using std::pow;

        return real_t(4./5) * K_0<real_t>() * T / p / (
          pow(
            real_t(2) * R_d<real_t>() * T
            / si::metres_per_second / si::metres_per_second, 
            real_t(.5)
          ) * si::metres_per_second
        );
      }
    };
  };
};
