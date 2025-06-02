#pragma once

#include "units.hpp"

namespace libcloudphxx
{
  namespace common
  {
    namespace ice_nucleation
    {
      enum class INP_t {mineral}; // types of ice nucleating particles, TODO: add more types

      // Freezing temperature as defined in Shima et al., 2020
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::temperature, real_t> T_freez(
      const INP_t& INP_type,                          // type of ice nucleating particle
      const quantity<si::length, real_t> rd_insol,    // radius of ice nucleating (insoluble) particle
      int rng_seed                                    // seed for random number generator
        ) {
        real_t A = real_t(4) * pi<real_t>() * rd_insol * rd_insol / si::square_meters; // surface area of the insoluble particle

        std::mt19937 gen(rng_seed); // random number generator
        std::uniform_real_distribution<> dis(0, 1);
        real_t Y = dis(gen); // random number between [0, 1]

        // Sampling T_freez from its CDF using inverse transform sampling
        if (INP_type == INP_t::mineral)
        {
          return real_t(273.15) + (real_t(8.934) - std::log(- std::log(1-Y) / A) ) / real_t(0.517);
        }
      }

    };
  };
};
