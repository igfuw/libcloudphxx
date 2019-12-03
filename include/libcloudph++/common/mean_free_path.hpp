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
#if !defined(__NVCC__)
        //using std::sqrt;
#endif
        return real_t(2) * (D_0<real_t>() / si::metres_per_second) / (
          real_t(sqrt(
            real_t(2) * real_t((R_v<real_t>() / si::joules * si::kilograms) * T) 
            ))
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
//using std::pow;

        return real_t(.8) * (K_0<real_t>() * T / p / si::metres_per_second) / (
          real_t(sqrt(
            real_t(2) * real_t((R_d<real_t>() / si::joules * si::kilograms) * T)
            ))
          );
      }
    };
  };
};
