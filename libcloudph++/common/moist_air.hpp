#pragma once

#include <libcloudph++/common/units.hpp>
#include <libcloudph++/common/macros.hpp>

namespace libcloudphxx
{
  namespace common
  {
    namespace moist_air
    {
      // specific heat capacities
      libcloudphxx_declare_const_macro(c_pd, 1005, si::joules / si::kilograms / si::kelvins) // dry air
      libcloudphxx_declare_const_macro(c_pv, 1850, si::joules / si::kilograms / si::kelvins) // water vapour
      libcloudphxx_declare_const_macro(c_pw, 4218, si::joules / si::kilograms / si::kelvins) // liquid water

      // molar masses
      libcloudphxx_declare_const_macro(M_d, 0.02896, si::kilograms / si::moles) // dry air
      libcloudphxx_declare_const_macro(M_v, 0.01802, si::kilograms / si::moles) // water vapour
      libcloudphxx_derived_const_macro(eps, M_v<real_t>() / M_d<real_t>()) // aka epsilon

      // universal gas constant (i.e. the Boltzmann times the Avogadro constants)
      libcloudphxx_declare_const_macro(kaBoNA, 8.314472, si::joules / si::kelvins / si::moles)

      // gas constants
      libcloudphxx_derived_const_macro(R_d, kaBoNA<real_t>() / M_d<real_t>()) // dry air
      libcloudphxx_derived_const_macro(R_v, kaBoNA<real_t>() / M_v<real_t>()) // water vapour

      // Exner function exponent for dry air
      libcloudphxx_derived_const_macro(R_d_over_c_pd, R_d<real_t>() / c_pd<real_t>())

      // water density
      libcloudphxx_declare_const_macro(rho_w, 1e3, si::kilograms / si::cubic_metres)

      // mixing rule for extensive quantitites (i.e. using mass mixing ratio)
      template <typename real_t, typename quant>
      auto constexpr mix(
	const quant &dry, 
	const quant &vap, 
	const quantity<si::dimensionless, real_t> &r
      ) libcloudphxx_decltype_return_macro(
	(dry + r * vap) / (1 + r)
      )

      // gas constant for moist air
      libcloudphxx_declare_funct_macro auto R(
	const quantity<si::dimensionless, real_t> &r
      ) libcloudphxx_decltype_return_macro(
	mix(R_d<real_t>(), R_v<real_t>(), r)
      )
     
      // specific heat capacity of moist air
      libcloudphxx_declare_funct_macro auto c_p(
	const quantity<si::dimensionless, real_t> &r
      ) libcloudphxx_decltype_return_macro(
	mix(c_pd<real_t>(), c_pv<real_t>(), r)
      )

      // water vapour partial pressure as a function of mixing ratio
      libcloudphxx_declare_funct_macro quantity<si::pressure, real_t> p_v(
	const quantity<si::pressure, real_t> &p,
	const quantity<si::dimensionless, real_t> &r
      )
      {
	return p * r / (r + eps<real_t>());
      }
    };
  };
};
