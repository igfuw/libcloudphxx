#pragma once

#include "const_cp.hpp"
#include "theta_std.hpp"

// theta dry: \theta = T * (p_1000 / p_dry)^{R_d / c_{pd}}
// theta std: \theta = T * (p_1000 / p    )^{R_d / c_{pd}}

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

      template <typename real_t>
      BOOST_GPU_ENABLED
      inline quantity<si::temperature, real_t> T(
        const quantity<si::temperature, real_t> &th, // theta dry!!!
        const quantity<si::mass_density, real_t> &rhod
      ) {

        /*
                return si::kelvins * pow(
          th / si::kelvins
          * pow(rhod * R_d<real_t>() / p_1000<real_t>() * si::kelvins, R_d<real_t>() / c_pd<real_t>()),
          c_pd<real_t>() / (c_pd<real_t>() - R_d<real_t>())
        );
                */
        return pow(
          real_t(th / si::kelvins)
          * pow(real_t(rhod * R_d<real_t>() / p_1000<real_t>() * si::kelvins), real_t(R_d<real_t>() / c_pd<real_t>())), 
          real_t(c_pd<real_t>() / (c_pd<real_t>() - R_d<real_t>()))
        ) * si::kelvins;
      }

      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::pressure, real_t> p(
        const quantity<si::mass_density, real_t> &rhod,
	const quantity<si::dimensionless, real_t> &r,
	const quantity<si::temperature, real_t> &T
      ) {
        return rhod         * (R_d<real_t>() + r * R_v<real_t>()) * T;
        //     ^^^^           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   
        //     rho/(1+r)      R(r)*(1+r)
      }

      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::temperature, real_t> d_th_d_rv(
	const quantity<si::temperature, real_t> &T,
	const quantity<si::temperature, real_t> &th // theta dry!!!
      ) {
	return - th / T * const_cp::l_v<real_t>(T) / c_pd<real_t>();
      }

      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::temperature, real_t> std2dry(
        const quantity<si::temperature, real_t> &th_std,
        const quantity<si::dimensionless, real_t> &r
      ) {
#if !defined(__NVCC__)
        //using std::pow;
#endif
        return th_std * pow(
          1 + r * R_v<real_t>() / R_d<real_t>(), 
          R_d<real_t>() / c_pd<real_t>()
        );
      }

      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::temperature, real_t> dry2std(
        const quantity<si::temperature, real_t> &th_dry,
        const quantity<si::dimensionless, real_t> &r
      ) {
#if !defined(__NVCC__)
        //using std::pow;
#endif
        return th_dry / pow(
          1 + r * R_v<real_t>() / R_d<real_t>(), 
          R_d<real_t>() / c_pd<real_t>()
        );
      }
    };
  };
};
