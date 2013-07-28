#pragma once

#include <libcloudph++/common/const_cp.hpp>
#include <libcloudph++/common/theta_std.hpp>

// theta dry: \theta = (p_1000 / p_dry)^{R_d / c_{pd}}
// theta std: \theta = (p_1000 / p)^{R_d / c_{pd}}

namespace libcloudphxx
{
  namespace common
  {
    namespace theta_dry
    {
      using moist_air::R;
      using moist_air::R_d;
      using moist_air::R_v;
      using moist_air::c_pd;
      using moist_air::p_v;
      using theta_std::p_1000;

      // TODO: eq. no
      template <typename real_t>
      BOOST_GPU_ENABLED
      inline quantity<si::temperature, real_t> T(
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
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::pressure, real_t> p(
        const quantity<si::mass_density, real_t> &rhod,
	const quantity<si::dimensionless, real_t> &r,
	const quantity<si::temperature, real_t> &T
      ) {
        return rhod * (R_d<real_t>() + r * R_v<real_t>()) * T;
      }
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::pressure, real_t> p(
        const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> &rhod_T,
	const quantity<si::dimensionless, real_t> &r
      ) {
        return rhod_T * (R_d<real_t>() + r * R_v<real_t>());
      }


      // see eq. TODO
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> d_rhodtheta_d_rv(
	const quantity<si::temperature, real_t> &T,
	const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> &rhod_th // theta dry!!!
      ) {
	return - rhod_th * const_cp::l_v<real_t>(T) / c_pd<real_t>() / T;
      }
    };
  };
};
