#pragma once

#include <libcloudph++/common/const_cp.hpp>

namespace libcloudphxx
{
  namespace common
  {
    namespace theta
    {
      // TODO: we use dry air here!
      using moist_air::R;
      using moist_air::R_d;
      using moist_air::eps;
      using moist_air::R_d_over_c_pd;
      using moist_air::c_p;
      using moist_air::p_v;

      // pressure in the definition of potential temperature
      libcloudphxx_declare_const_macro(p_1000, 100000, si::pascals)

      // Exner function exponent for moist air
      libcloudphxx_declare_funct_macro quantity<si::dimensionless, real_t> R_over_c_p(
	const quantity<si::dimensionless, real_t> &r
      )
      {
	return R<real_t>(r) / c_p<real_t>(r);
      }


      // Exner function for dry air
      libcloudphxx_declare_funct_macro quantity<si::dimensionless, real_t> exner(
	const quantity<si::pressure, real_t> &p
      )
      {
	return pow(p / p_1000<real_t>(), R_d_over_c_pd<real_t>());
      }


      // Exner function for moist air
      libcloudphxx_declare_funct_macro quantity<si::dimensionless, real_t> exner(
	const quantity<si::pressure, real_t> &p,
	const quantity<si::dimensionless, real_t> &r
      )
      {
	return pow(p / p_1000<real_t>(), R_over_c_p<real_t>(r));
      }


      // dry air density as a function of p, theta and rv
      libcloudphxx_declare_funct_macro quantity<si::mass_density, real_t> rhod(
	const quantity<si::pressure, real_t> &p,
	const quantity<si::temperature, real_t> &th, // 
	const quantity<si::dimensionless, real_t> &rv
      )
      {
	return (p - p_v<real_t>(p, rv)) /
	  (exner<real_t>(p, rv) * R_d<real_t>() * th);
      }


      // temperature as a function theta, pressure and water vapour mixing ratio for moist air
      libcloudphxx_declare_funct_macro quantity<si::temperature, real_t> T(
	const quantity<si::temperature, real_t> &th, 
	const quantity<si::pressure, real_t> &p, 
	const quantity<si::dimensionless, real_t> &r 
      )
      {
	return th * exner<real_t>(p, r);
      }


      // pressure as a function of "theta times dry air density" and water vapour mixing ratio
      libcloudphxx_declare_funct_macro quantity<si::pressure, real_t> p(
	const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> rhod_th,
	const quantity<si::dimensionless, real_t> &r 
      )
      {
	return p_1000<real_t>() * real_t(pow(
	  (rhod_th * R_d<real_t>())
	    / p_1000<real_t>() * (real_t(1) + r / eps<real_t>()),
	  1 / (1 - R_over_c_p(r))
	));
      }

      // dtheta^star_drv from First Law for theta^star
      libcloudphxx_declare_funct_macro quantity<si::temperature, real_t> dtheta_drv(
	const quantity<si::temperature, real_t> &T,
	const quantity<si::pressure, real_t> &p,
	const quantity<si::dimensionless, real_t> &r,
	const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> &rhod_th,
	const quantity<si::mass_density, real_t> &rhod
      )
      {
	return - rhod_th / rhod * (
	  // the 'liquid water' term
	  const_cp::l_v<real_t>(T)
	    / real_t(pow(1 + r, 2))
	    / c_p(r)
	    / T
	);
      }
    };
  };
};
