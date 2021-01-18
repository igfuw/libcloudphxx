#pragma once

#include "units.hpp" // TODO: do detail?
#include "macros.hpp" // TODO: do detail?
#include "moist_air.hpp"

namespace libcloudphxx
{
  namespace common
  {
    namespace const_cp
    {
      using moist_air::c_pw;
      using moist_air::c_pv;
      using moist_air::R_v;
      using moist_air::eps;
      typedef divide_typeof_helper<si::energy, si::mass>::type energy_over_mass;


      // water triple point parameters
      libcloudphxx_const(si::pressure, p_tri, 611.73, si::pascals) // pressure
      libcloudphxx_const(si::temperature, T_tri, 273.16, si::kelvins) // temperature
      libcloudphxx_const(energy_over_mass, l_tri, 2.5e6, si::joules / si::kilograms) // latent heat of evaporation

      // saturation vapour pressure for water assuming constant c_p_v and c_p_w
      // with constants taken at triple point
      // (solution to the Clausius-Clapeyron equation assuming rho_vapour << rho_liquid)
//<listing-1>
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::pressure, real_t> p_vs(
        const quantity<si::temperature, real_t> &T
      ) 
//</listing-1>
      {
        return p_tri<real_t>() * exp(
          (l_tri<real_t>() + (c_pw<real_t>() - c_pv<real_t>()) * T_tri<real_t>()) / R_v<real_t>() * (real_t(1) / T_tri<real_t>() - real_t(1) / T)
          - (c_pw<real_t>() - c_pv<real_t>()) / R_v<real_t>() * std::log(T / T_tri<real_t>())
        );
      }

      // saturation vapour mixing ratio for water as a function of pressure and temperature
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::dimensionless, real_t> r_vs(
        const quantity<si::temperature, real_t> &T,
        const quantity<si::pressure, real_t> &p
      ) {
        return eps<real_t>() / (p / p_vs<real_t>(T) - 1);
      }

      // latent heat for constant c_p
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<divide_typeof_helper<si::energy, si::mass>::type , real_t> l_v(
        const quantity<si::temperature, real_t> &T
      ) {
        return l_tri<real_t>() + (c_pv<real_t>() - c_pw<real_t>()) * (T - T_tri<real_t>());
      }
    };
  };
};
