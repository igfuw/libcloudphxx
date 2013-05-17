#pragma once

#include <libcloudph++/common/const_cp.hpp>

// theta dry: \theta = (p_1000 / p_dry)^{R_d / c_{pd}}
// normal theta : \theta = (p_1000 / p)^{R_d / c_{pd}}

namespace libcloudphxx
{
  namespace common
  {
    namespace theta
    {
      using moist_air::R;
      using moist_air::R_d;
      using moist_air::R_v;
      using moist_air::R_d_over_c_pd;
      using moist_air::c_pd;
      using moist_air::p_v;

      // pressure in the definition of potential temperature
      libcloudphxx_declare_const_macro(p_1000, 100000, si::pascals)

      // dry air density as a function of p, theta and rv
      libcloudphxx_declare_funct_macro quantity<si::mass_density, real_t> rhod(
	const quantity<si::pressure, real_t> &p,
	const quantity<si::temperature, real_t> &th, // normal theta !!! TODO: document / separate file?
	const quantity<si::dimensionless, real_t> &rv
      )
      {
	return (p - p_v<real_t>(p, rv)) /
	  (pow(p / p_1000<real_t>(), R_d_over_c_pd<real_t>()) * R_d<real_t>() * th);
      }


      // TODO: eq. no
      libcloudphxx_declare_funct_macro quantity<si::temperature, real_t> T(
        const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> &rhod_th, // theta dry!!!
        const quantity<si::mass_density, real_t> &rhod
      ) {
        return si::kelvins * pow(
          rhod_th / rhod / si::kelvins
          * pow(rhod * R_d<real_t>() / p_1000<real_t>() * si::kelvins, R_d<real_t>() / c_pd<real_t>()), 
          c_pd<real_t>() / (c_pd<real_t>() - R_d<real_t>())
        );
      }


      // TODO: eq. 
      libcloudphxx_declare_funct_macro quantity<si::pressure, real_t> p(
        const quantity<si::mass_density, real_t> &rhod,
	const quantity<si::dimensionless, real_t> &r,
	const quantity<si::temperature, real_t> &T
      )
      {
        return rhod * (R_d<real_t>() + r * R_v<real_t>()) * T;
      }


      // see eq. TODO
      libcloudphxx_declare_funct_macro quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> d_rhodtheta_d_rv(
	const quantity<si::temperature, real_t> &T,
	const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> &rhod_th // theta dry!!!
      )
      {
	return - rhod_th * const_cp::l_v<real_t>(T) / c_pd<real_t>() / T;
      }
    };
  };
};
