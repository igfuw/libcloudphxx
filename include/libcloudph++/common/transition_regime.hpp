#pragma once

#include "units.hpp"
#include "macros.hpp"

namespace libcloudphxx
{
  namespace common
  {
    namespace transition_regime
    {
      // see Laaksonen et al. 2005 (after Fuchs and Sutugin)
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::dimensionless, real_t> beta(
        quantity<si::dimensionless, real_t> Kn // Knudsen number
      ) 
      {
        return (1 + Kn) / (1 + real_t(1.71) * Kn + 1.33 * Kn*Kn);
      }
    };
  };
};
