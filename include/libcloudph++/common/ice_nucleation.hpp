#pragma once

#include "units.hpp"
#include "const_cp.hpp"

namespace libcloudphxx
{
  namespace common
  {
    namespace ice_nucleation
    {
      enum class INP_t {mineral}; // types of ice nucleating particles, TODO: add more types

      // Inverse CDF for singular freezing temperature as defined in eq. 1 in Shima et al., 2020
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::temperature, real_t> T_freeze_CDF_inv(
      const INP_t& INP_type,      // type of ice nucleating particle
      const real_t rd2_insol,     // radius squared of insoluble particle in m^2
      const real_t rand           // random number between [0, 1]
        ) {
        real_t A = real_t(4)
        #if !defined(__NVCC__)
            * pi<real_t>()
        #else
            * CUDART_PI
        #endif
        * rd2_insol; // surface area of the insoluble particle

        if (INP_type == INP_t::mineral && A > real_t(1e-20))
        {
          return (real_t(273.15) + (real_t(8.934) - log(- log(real_t(1.) - rand) / A) ) / real_t(0.517)) * si::kelvin;
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
          const real_t &rd2_insol = thrust::get<0>(tpl);  // from rd2 vector
          const real_t &rand         = thrust::get<1>(tpl);  // from rand vector

          return ice_nucleation::template T_freeze_CDF_inv<real_t>(
            INP_type,
            rd2_insol,
            rand
          ).value();
        }
      };

      // Probability of time-dependent freezing as in Arabas et al., 2025
      template <typename real_t>
      BOOST_GPU_ENABLED
      real_t p_freeze(
      const INP_t& INP_type,     // type of ice nucleating particle
      const real_t rd2_insol,    // radius squared of insoluble particle in m^2
      const real_t T,            // temperature in kelvin
      const real_t dt            // time step in seconds
        )
      {
        if (rd2_insol > real_t(0))
        {
          real_t A = real_t(4)
          #if !defined(__NVCC__)
              * pi<real_t>()
          #else
              * CUDART_PI
          #endif
          * rd2_insol; // surface area of the insoluble particle
          real_t d_aw = real_t(1) - const_cp::p_vsi<real_t>(T * si::kelvin)/ const_cp::p_vs<real_t>(T * si::kelvin); // water activity
          if (INP_type == INP_t::mineral)
          {
            real_t J = pow(real_t(10), real_t(-1.35) + real_t(22.62) * d_aw) * real_t(1e4); // nucleation rate
            return 1 - exp(- J * A * dt);
          }
          else
            return real_t(0.); // TODO: other INP types
        }
        else
          return T > real_t(235.15) ? real_t(0) : real_t(1); // homogeneous freezing at -38 C
      }


      template<typename real_t>
      struct p_freeze_functor
      {
        INP_t INP_type;
        real_t dt;

        p_freeze_functor(INP_t INP_type, real_t dt)
          : INP_type(INP_type), dt(dt) {}

        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t> &tpl) const
        {
          const real_t &rd2_insol = thrust::get<0>(tpl);  // radius squared of insoluble particle
          const real_t &T         = thrust::get<1>(tpl);  // temperature in kelvin

          return ice_nucleation::p_freeze<real_t>(
            INP_type,
            rd2_insol,
            T,
            dt
          );
        }
      };

    };
  };
};
