#pragma once

#include "units.hpp"

namespace libcloudphxx
{
  namespace common
  {
    namespace ice_nucleation
    {
      enum class INP_t {mineral}; // types of ice nucleating particles, TODO: add more types

      // Inverse CDF for freezing temperature as defined in eq. 1 in Shima et al., 2020
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::temperature, real_t> T_freeze_CDF_inv(
      const INP_t& INP_type,                          // type of ice nucleating particle
      const quantity<si::length, real_t> rd3_insol,    // radius cubed of ice nucleating (insoluble) particle
      const real_t rand                                    // random number between [0, 1]
        ) {
        real_t A = real_t(4) * pi<real_t>() * std::pow(rd3_insol/si::meters, 2/3); // surface area of the insoluble particle

        if (INP_type == INP_t::mineral && A > std::numeric_limits<real_t>::epsilon())
        {
          return (real_t(273.15) + (real_t(8.934) - std::log(- std::log(1 - rand) / A) ) / real_t(0.517)) * si::kelvins;
        }
        else
        {
          return real_t(235.15) * si::kelvin; // the default freezing temperature is -38 C
        }
      }


      template<typename real_t>
      struct T_freeze_CDF_inv_functor
      {
        INP_t INP_type;

        T_freeze_CDF_inv_functor(INP_t INP_type)
          : INP_type(INP_type) {}

        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t> &tpl) const
        {
          const real_t &rd3_insol = thrust::get<0>(tpl);  // from rd3 vector
          const real_t &rand         = thrust::get<1>(tpl);  // from rand vector

          return ice_nucleation::T_freeze_CDF_inv<real_t>(
            INP_type,
            rd3_insol * si::meters,
            rand
          ).value();
        }
      };

    };
  };
};
