#pragma once

#include "units.hpp"
#include "const_cp.hpp"
#include <boost/math/special_functions/asin.hpp>

#if defined(__NVCC__)
#  include <math_constants.h>
#endif

namespace libcloudphxx
{
  namespace common
  {
    namespace ice_deposition
    {

      // capacitance of ice crystals
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::dimensionless, real_t> ice_capacitance(
      const quantity<si::length, real_t> ice_a,         // equatorial radius
      const quantity<si::length, real_t> ice_c          // polar radius
      )
      {
        if (ice_a < std::numeric_limits::epsilon) throw std::runtime_error("ice_a too small for calculation of capacitance");
        real_t phi = ice_c / ice_a;
        real_t e = std::sqrt(real_t(1) - phi * phi); // eccentricity
        return phi < real_t(1) ?
          ice_a * e / std::asin(e) : ice_a * e / std::log(1 + e) / phi;
      }


      // rate of change of ice mass
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<divide_typeof_helper<si::area, si::time>::type, real_t> dmi_dt(
        const quantity<diffusivity, real_t> D,            // D
        const quantity<thermal_conductivity, real_t> K,   // K
        const quantity<si::mass_density, real_t> rho_v,   // ambient water vapour density
        const quantity<si::temperature, real_t> T,        // ambient temperature
        const quantity<si::pressure, real_t> p,           // ambient pressure
        const quantity<si::dimensionless, real_t> RH_i,   // p_v/p_vsi = relative humidity w.r.t. ice
        const quantity<si::length, real_t> ice_a,         // ice equatorial radius
        const quantity<si::length, real_t> ice_c          // ice polar radius
      )
      {
        using moist_air::R_v;

        quantity<divide_typeof_helper<si::energy, si::mass>::type, real_t>
          l_s = const_cp::l_s<real_t>(T);

        return real_t(4) * pi<real_t>() * ice_capacitance(ice_a, ice_c)
          * (real_t(1) - real_t(1) / RH_i)
          / (
            real_t(1)
              / D
              / rho_v
            +
            l_s
              / K
              / RH_i
              / T
              * (l_s / R_v<real_t>() / T - real_t(1))
          )
        ;
      }


      // deposition density from Shima et al. (2020)
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::mass_density, real_t> rho_dep(
        const quantity<diffusivity, real_t> D,            // D
        const quantity<thermal_conductivity, real_t> K,   // K
        const quantity<si::mass_density, real_t> rho_v,   // ambient water vapour density
        const quantity<si::temperature, real_t> T,        // ambient temperature
        const quantity<si::dimensionless, real_t> RH_i,   // p_v/p_vsi = relative humidity w.r.t. ice
        const quantity<si::length, real_t> ice_a         // ice equatorial radius
      )
      {

        return (ice_growth_ratio(T) < real_t(1) && ice_a < real_t(1e-4)*si::meters) ?
              moist_air::rho_i<>() : rho_dep_CL94(D, K, rho_v, T, RH_i);
      }


      // deposition density from Chen and Lamb (1994)
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::mass_density, real_t> rho_dep_CL94(
        const quantity<diffusivity, real_t> D,            // D
        const quantity<thermal_conductivity, real_t> K,   // K
        const quantity<si::mass_density, real_t> rho_v,   // ambient water vapour density
        const quantity<si::temperature, real_t> T,        // ambient temperature
        const quantity<si::dimensionless, real_t> RH_i   // p_v/p_vsi = relative humidity w.r.t. ice
      )
      {
        using moist_air::R_v;
        quantity<divide_typeof_helper<si::energy, si::mass>::type, real_t>
          l_s = const_cp::l_s<real_t>(T);

        using const_cp::p_vs, const_cp::p_vsi;

        quantity<si::mass_density, real_t> d_rhoi = (std::min(RH_i, p_vs/p_vsi) - real_t(1))
          / RH_i
          / D
          / (
            real_t(1)
              / D
              / rho_v
            +
            l_s
              / K
              / RH_i
              / T
              * (l_s / R_v<real_t>() / T - real_t(1))
          );

        return moist_air::rho_i<>()
        * std::exp(- real_t(3) * std::max(d_rhoi/si::kilograms*si::cubic_meters - real_t(5e-5), real_t(0))
          / ice_growth_ratio(T));
      }


      // inherent ice growth ratio
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::dimensionless, real_t> ice_growth_ratio(const quantity<si::temperature, real_t> T)
      {
      }


    };
  };
};
