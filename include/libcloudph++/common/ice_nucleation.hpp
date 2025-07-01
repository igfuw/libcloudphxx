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
      quantity<si::temperature, real_t> T_freeze(
      const INP_t& INP_type,                          // type of ice nucleating particle
      const quantity<si::length, real_t> rd3_insol,    // radius cubed of ice nucleating (insoluble) particle
      int rng_seed                                    // seed for random number generator
        ) {
        real_t A = real_t(4) * pi<real_t>() * std::pow(rd3_insol/si::meters, 2/3); // surface area of the insoluble particle

        std::mt19937 gen(rng_seed); // random number generator
        std::uniform_real_distribution<> dis(0, 1);
        real_t Y = dis(gen); // random number between [0, 1]

        // Sampling T_freez from its CDF using inverse transform sampling
        if (A > std::numeric_limits<real_t>::epsilon && INP_type == INP_t::mineral)
        {
          return (real_t(273.15) + (real_t(8.934) - std::log(- std::log(1-Y) / A) ) / real_t(0.517)) * si::kelvins;
        }
        else
        {
          return real_t(235.15) * si::kelvin; // if rd_insol = 0  or INP type is unknown, the default freezing temperature is -38 C
        }
      }


      template<typename real_t>
      struct T_freeze_functor
      {
        INP_t INP_type;
        int rng_seed;
        T_freeze_functor(INP_t INP_type, int rng_seed)
          : INP_type(INP_type), rng_seed(rng_seed) {}
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &rd3_insol) const
        {
          return ice_nucleation::T_freeze<real_t>(INP_type, rd3_insol * si::meters, rng_seed).value();
        }
      };

    };
  };
};
