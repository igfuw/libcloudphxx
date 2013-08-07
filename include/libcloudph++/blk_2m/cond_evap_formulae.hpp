/** @file
  * @copyright University of Warsaw
  * @brief double-moment bulk condensation/evaporation parameterisation formulae
  * @date August 2013
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once
#include <libcloudph++/common/moist_air.hpp>

namespace libcloudphxx
{
  namespace blk_2m
  {
    namespace formulae
    {

      using namespace common::moist_air;
      using namespace common::const_cp;

      //vapour diffusivity in air (see Properties of air, Tracy, Welch, Porter)
      typedef divide_typeof_helper<si::area, si::time>::type diffusion;
      template<typename real_t>
      quantity<diffusion, real_t> air_vap_diff(
        quantity<si::temperature, real_t> T, 
        quantity<si::pressure, real_t> pres,
        quantity<diffusion, real_t> D_0 = 2.26 * 1e-5 *si::square_metres/si::seconds,
        quantity<si::pressure, real_t> p_0 = 100000 * si::pascal,
        quantity<si::temperature, real_t> T_0 = 273.15 * si::kelvin 
      ) {
        return D_0 * pow(T / T_0, 1.81) * (p_0 / pres);
      }

      // relaxation time for condensation/evaporation term for cloud and rain droplets
      // (see Khvorostyaov at al 2001 eq. 5)
      template<typename real_t>
      quantity<si::time, real_t> tau_relax(
        quantity<si::temperature, real_t> T, 
        quantity<si::pressure, real_t> p,
        quantity<si::length, real_t> r, //droplet radius
        quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> N
      ) {
        return 1. / (4 * pi<real_t>() * air_vap_diff(T, p) * N * r);
      }
 
      //drv_s/dT (derived from Clapeyron equation and pv = rv * rho_d * R_v * T)
      typedef divide_typeof_helper<si::dimensionless, si::temperature>::type one_over_temperature;
      template<typename real_t> 
      quantity<one_over_temperature, real_t> drv_s_dT(
        quantity<si::temperature, real_t> T,
        quantity<si::pressure, real_t> p
      ) {
        return l_v(T) * r_vs(T, p) / R_v<real_t>() / pow<2>(T);
      }

      //condensation/evaporation rate
      template<typename real_t>
      quantity<si::frequency, real_t> cond_evap_rate(
        quantity<si::temperature, real_t> T, 
        quantity<si::pressure, real_t> p,
        quantity<si::length, real_t> r, //droplet radius
        quantity<si::dimensionless, real_t> r_v,
        quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> N
      ) {
        return (r_v-r_vs(T,p)) / tau_relax(T, p, r, N) / (1 + drv_s_dT(T, p) * l_v(T) / c_p(r_v));
      }                                                                     //TODO check ^ 

    };
  };
};
