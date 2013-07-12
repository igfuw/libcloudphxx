#pragma once

#include <libcloudph++/common/units.hpp> // TODO: do detail?
#include <libcloudph++/common/macros.hpp> // TODO: do detail?
#include <libcloudph++/common/moist_air.hpp>

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

      // water triple point parameters
      libcloudphxx_declare_const_macro(p_tri, 611.73, si::pascals) // pressure
      libcloudphxx_declare_const_macro(T_tri, 273.16, si::kelvins) // temperature
      libcloudphxx_declare_const_macro(l_tri, 2.5e6, si::joules / si::kilograms) // latent heat of evaporation

      // saturation vapour pressure for water assuming constant c_p_v and c_p_w
      // with constants taken at triple point
      // (solution to the Clausius-Clapeyron equation assuming rho_vapour << rho_liquid)
      libcloudphxx_declare_funct_macro quantity<si::pressure, real_t> p_vs(
	const quantity<si::temperature, real_t> &T
      )
      {
	return p_tri<real_t>() * exp(
	  (l_tri<real_t>() + (c_pw<real_t>() - c_pv<real_t>()) * T_tri<real_t>()) / R_v<real_t>() * (real_t(1) / T_tri<real_t>() - real_t(1) / T)
	  - (c_pw<real_t>() - c_pv<real_t>()) / R_v<real_t>() * log(T / T_tri<real_t>())
	);
      }

      // saturation vapour mixing ratio for water as a function of pressure and temperature
      libcloudphxx_declare_funct_macro quantity<si::dimensionless, real_t> r_vs(
	const quantity<si::temperature, real_t> &T,
	const quantity<si::pressure, real_t> &p
      )
      {
	return eps<real_t>() / (p / p_vs<real_t>(T) - 1);
      }

      // latent heat for constant c_p
      libcloudphxx_declare_funct_macro quantity<divide_typeof_helper<si::energy, si::mass>::type , real_t> l_v(
	const quantity<si::temperature, real_t> &T
      )
      {
	return l_tri<real_t>() + (c_pv<real_t>() - c_pw<real_t>()) * (T - T_tri<real_t>());
      }
    };
  };
};
