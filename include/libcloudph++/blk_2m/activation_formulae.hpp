/** @file
  * @copyright University of Warsaw
  * @brief double-moment bulk activation parameterisation formulae
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libcloudph++/common/moist_air.hpp>
#include <libcloudph++/common/kelvin_term.hpp>
#include <libcloudph++/common/const_cp.hpp>
#include <libcloudph++/common/earth.hpp>

namespace libcloudphxx
{
  namespace blk_2m
  {
    namespace formulae
    {
      using common::earth::rho_stp;

      //beta parameter for defining solubility of dry aerosol (see Khvorostyanov and Curry 2006)
      //beta=0.5 means that soluble fraction is proportional to the volume of the ary aerosol
      libcloudphxx_const(si::dimensionless, beta_default, .5, 1)

      //in activation term all activated droples are assumed to have the radius of 1 um
      libcloudphxx_const_derived(si::mass, ccnmass, 
        real_t(4./3 * pi<real_t>()) * pow<3>(1e-6 * si::metres) * common::moist_air::rho_w<real_t>());

      //helper for activation formulae (see eq. 11 in Morrison and Grabowski 2007)
      template <typename real_t>
      quantity<si::dimensionless, real_t> s_0(
        quantity<si::temperature,   real_t> T,
        quantity<si::length,        real_t> mean_rd,
        quantity<si::dimensionless, real_t> chem_b,
        quantity<si::dimensionless, real_t> beta = beta_default<real_t>()
      ) {
        return pow(mean_rd / si::metres, -(1+beta)) * sqrt(4 * pow(common::kelvin::A<real_t>(T) / si::metres, 3) / 27 / chem_b);
      }           // can't have pow<3/2>(si:: )

      //helper for activation formulae (see eq. 10 in Morrison and Grabowski 2007)
      template <typename real_t>
      quantity<si::dimensionless, real_t> s(
        quantity<si::pressure,     real_t> p,
        quantity<si::temperature,  real_t> T, 
        quantity<si::mass_density, real_t> rhod,
        quantity<si::mass_density, real_t> rhod_rv
      ) {
        return rhod_rv / rhod / common::const_cp::r_vs<real_t>(T, p) - real_t(1);
      }

      //helper for activation formulae (see eq. 12 in Morrison and Grabowski 2007)
      template<typename real_t>
      quantity<si::dimensionless, real_t> sdev_rd_s(
        quantity<si::dimensionless, real_t> sdev_rd,
        quantity<si::dimensionless, real_t> beta = beta_default<real_t>()
      ) {
        return pow(sdev_rd, 1+beta);
      }

      //helper for activation formulae (see eq. 10 in Morrison and Grabowski 2007)
      template <typename real_t>
      quantity<si::dimensionless, real_t> u(
        quantity<si::pressure,      real_t> p,
        quantity<si::temperature,   real_t> T,
        quantity<si::mass_density,  real_t> rhod,
        quantity<si::mass_density,  real_t> rhod_rv,
        quantity<si::length,        real_t> mean_rd,
        quantity<si::dimensionless, real_t> sdev_rd,
        quantity<si::dimensionless, real_t> chem_b,
        quantity<si::dimensionless, real_t> beta = beta_default<real_t>()
      ) {
        return log(s_0(T, mean_rd, chem_b) / s(p, T, rhod, rhod_rv)) / sqrt(2) / log(sdev_rd_s(sdev_rd));
      }

      typedef divide_typeof_helper<si::dimensionless, si::volume>::type one_over_volume;
      //helper for activation formulae (see eq. 10 in Morrison and Grabowski 2007)
      template <typename real_t>
      quantity<one_over_volume, real_t> N_c_p(
        quantity<si::pressure,      real_t> p,
        quantity<si::temperature,   real_t> T,
        quantity<si::mass_density,  real_t> rhod,
        quantity<si::mass_density,  real_t> rhod_rv,
        quantity<si::length,        real_t> mean_rd,
        quantity<si::dimensionless, real_t> sdev_rd,
        quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> N_stp,
        quantity<si::dimensionless, real_t> chem_b,
        quantity<si::dimensionless, real_t> beta = beta_default<real_t>()
      ) {
        return (N_stp / rho_stp<real_t>() * rhod) / real_t(2.) * erfc(u(p, T, rhod, rhod_rv, mean_rd, sdev_rd, chem_b)); 
      }

      typedef divide_typeof_helper<one_over_volume, si::time>::type one_over_volume_time;
      //activation formulae (see eq. 13 in Morrison and Grabowski 2007)
      template <typename real_t>
      quantity<one_over_volume_time, real_t> activation_rate(
        quantity<si::pressure,      real_t> p,
        quantity<si::temperature,   real_t> T,
        quantity<si::mass_density,  real_t> rhod,
        quantity<si::mass_density,  real_t> rhod_rv,
        quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> rhod_nc,
        quantity<si::length,        real_t> mean_rd,
        quantity<si::dimensionless, real_t> sdev_rd,
        quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> N_stp,
        const quantity<si::time, real_t> dt,
        quantity<si::dimensionless, real_t> chem_b,
        quantity<si::dimensionless, real_t> beta = beta_default<real_t>()
      ) {
        return std::max(
          real_t(0) / si::cubic_metres / si::seconds,
          (N_c_p<real_t>(p, T, rhod, rhod_rv, mean_rd, sdev_rd, N_stp, chem_b) - rhod_nc) / dt
        );
      }
 
    };
  };
};
