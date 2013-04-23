#pragma once
#include "phc.hpp"

namespace phc
{
  // saturation vapour pressure for water assuming constant c_p_v and c_p_w
  // with constants taken at triple point
  // (solution to the Clausius-Clapeyron equation assuming rho_vapour << rho_liquid)
  phc_declare_funct_macro quantity<si::pressure, real_t> p_vs(
    quantity<si::temperature, real_t> T
  )
  {
    return p_tri<real_t>() * exp(
      (l_tri<real_t>() + (c_pw<real_t>() - c_pv<real_t>()) * T_tri<real_t>()) / R_v<real_t>() * (real_t(1) / T_tri<real_t>() - real_t(1) / T)
      - (c_pw<real_t>() - c_pv<real_t>()) / R_v<real_t>() * log(T / T_tri<real_t>())
    );
  }

  // saturation vapour mixing ratio for water as a function of pressure and temperature
  phc_declare_funct_macro quantity<mixing_ratio, real_t> r_vs(
    quantity<si::temperature, real_t> T,
    quantity<si::pressure, real_t> p
  )
  {
    return eps<real_t>() / (p / p_vs<real_t>(T) - 1);
  }

  // latent heat for constant c_p
  phc_declare_funct_macro quantity<divide_typeof_helper<si::energy, si::mass>::type , real_t> l_v(
    quantity<si::temperature, real_t> T
  )
  {
    return l_tri<real_t>() + (c_pv<real_t>() - c_pw<real_t>()) * (T - T_tri<real_t>());
  }
};
