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
        real_t(4./3) * pi<real_t>() * std::pow(real_t(1e-6), real_t(3)) * si::cubic_metres * common::moist_air::rho_w<real_t>());

      //helper for activation formulae (see eq. 11 in Morrison and Grabowski 2007)
      template <typename real_t>
      inline quantity<si::dimensionless, real_t> s_0(
        const quantity<si::temperature,   real_t> &T,
        const quantity<si::length,        real_t> &mean_rd,
        const quantity<si::dimensionless, real_t> &chem_b,
        const quantity<si::dimensionless, real_t> beta = beta_default<real_t>()
      ) {
        return std::pow(mean_rd / si::metres, -(1+beta)) * std::sqrt(4 * std::pow(common::kelvin::A<real_t>(T) / si::metres, 3) / 27 / chem_b);
      }           // can't have std::pow<3/2>(si:: )

      // helper for activation formulae (see eq. 10 in Morrison and Grabowski 2007)
      template <typename real_t>
      inline quantity<si::dimensionless,  real_t> s(
        const quantity<si::pressure,      real_t> &p,
        const quantity<si::temperature,   real_t> &T, 
        const quantity<si::dimensionless, real_t> &rv
      ) {
        return rv / common::const_cp::r_vs<real_t>(T, p) - real_t(1);
      }

      // helper for activation formulae (see eq. 12 in Morrison and Grabowski 2007)
      template<typename real_t>
      inline quantity<si::dimensionless, real_t> sdev_rd_s(
        const quantity<si::dimensionless, real_t> &sdev_rd,
        const quantity<si::dimensionless, real_t> beta = beta_default<real_t>()
      ) {
        return std::pow(sdev_rd, 1+beta);
      }

      //helper for activation formulae (see eq. 10 in Morrison and Grabowski 2007)
      template <typename real_t>
      inline quantity<si::dimensionless, real_t> u(
        const quantity<si::pressure,      real_t> &p,
        const quantity<si::temperature,   real_t> &T,
        const quantity<si::dimensionless, real_t> &rv,
        const quantity<si::length,        real_t> &mean_rd,
        const quantity<si::dimensionless, real_t> &sdev_rd,
        const quantity<si::dimensionless, real_t> &chem_b,
        const quantity<si::dimensionless, real_t> &RH_max,
        const quantity<si::dimensionless, real_t> beta = beta_default<real_t>()
      ) {
        return log(
          s_0(T, mean_rd, chem_b) / 
          std::min(real_t(s(p, T, rv)), real_t(RH_max - 1))
        ) / std::sqrt(2) / log(sdev_rd_s(sdev_rd));
      }

      // helper for activation formulae (see eq. 10 in Morrison and Grabowski 2007)
      template <typename real_t>
      inline quantity<divide_typeof_helper<si::dimensionless, si::mass>::type, real_t> n_c_p(
        const quantity<si::pressure,      real_t> &p,
        const quantity<si::temperature,   real_t> &T,
        const quantity<si::dimensionless, real_t> &rv,
        const quantity<si::length,        real_t> &mean_rd,
        const quantity<si::dimensionless, real_t> &sdev_rd,
        const quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> &N_stp,
        const quantity<si::dimensionless, real_t> &chem_b,
        const quantity<si::dimensionless, real_t> &RH_max,
        const quantity<si::dimensionless, real_t> beta = beta_default<real_t>()
      ) {
        return (N_stp / rho_stp<real_t>()) / real_t(2.) * std::erfc(u(p, T, rv, mean_rd, sdev_rd, chem_b, RH_max)); 
      }

      // activation formulae (see eq. 13 in Morrison and Grabowski 2007)
      template <typename real_t>
      inline quantity<divide_typeof_helper<si::frequency, si::mass>::type, real_t> activation_rate(
        const quantity<divide_typeof_helper<si::dimensionless, si::mass>::type, real_t> &n_ccn,
        const quantity<divide_typeof_helper<si::dimensionless, si::mass>::type, real_t> &nc,
        const quantity<si::time, real_t> &dt
      ) {
        return std::max(
          real_t(0) / si::kilograms / si::seconds,
          (n_ccn - nc) / dt
        );
      }
 
    };
  };
};
