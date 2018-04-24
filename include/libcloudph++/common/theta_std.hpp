#pragma once

#include <libcloudph++/common/const_cp.hpp>

// theta dry: \theta = (p_1000 / p_dry)^{R_d / c_{pd}}
// theta std: \theta = (p_1000 / p)^{R_d / c_{pd}}

namespace libcloudphxx
{
  namespace common
  {
    namespace theta_std
    {
      using moist_air::R;
      using moist_air::R_d;
      using moist_air::c_pd;
      using moist_air::p_v;

      // pressure in the definition of potential temperature
      libcloudphxx_const(si::pressure, p_1000, 100000, si::pascals)

      // dry air density as a function of p, theta and rv
      template <typename real_t>
      quantity<si::mass_density, real_t> rhod(
	const quantity<si::pressure, real_t> &p,
	const quantity<si::temperature, real_t> &th_std, 
	const quantity<si::dimensionless, real_t> &rv
      ) {
	return (p - p_v<real_t>(p, rv)) /
	  (pow(p / p_1000<real_t>(), R_d<real_t>()/c_pd<real_t>()) * R_d<real_t>() * th_std);
      }

      // Exner pressure
      template <typename real_t>
      quantity<si::dimensionless, real_t> exner(
	const quantity<si::pressure, real_t> &p
      ) {
	return pow(p / p_1000<real_t>(), R_d<real_t>()/c_pd<real_t>());
      }
    };
  };
};
