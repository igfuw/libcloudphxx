/** @file
  * @copyright University of Warsaw
  * @brief double-moment bulk condensation/evaporation parameterisation formulae
  * @date August 2013
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once
#include <libcloudph++/common/moist_air.hpp>
#include <libcloudph++/common/const_cp.hpp>
#include <libcloudph++/common/ventil.hpp>
#include <libcloudph++/common/vterm.hpp>
#include <libcloudph++/common/earth.hpp>
#include <libcloudph++/blk_2m/terminal_vel_formulae.hpp>

namespace libcloudphxx
{
  namespace blk_2m
  {
    namespace formulae
    {
      using namespace common::moist_air;
      using namespace common::const_cp;
      using namespace common::ventil; //Schmidt number
      using namespace common::vterm;  //air viscosity
      using namespace common::earth;  //rho_stp

      // relaxation time for condensation/evaporation term for cloud droplets
      // (see eq.4 in Morrison 2005 but with f1=1 and f2=0 - neglecting ventilation coeffs for cloud droplets)
      // (or Khvorostyaov at al 2001 eq. 5)
      template<typename real_t>
      quantity<si::time, real_t> tau_relax_c(
        const quantity<si::temperature, real_t> T, 
        const quantity<si::pressure, real_t> p,
        const quantity<si::length, real_t> r, //droplet radius
        const quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> N
      ) {
        return real_t(1.) / (4 * pi<real_t>() * D_0<real_t>() * N * r);
      }

      //ventilation coefficients TODO - check are those really those numbers?
      //see Morrison 2005 eq.4
      //but also Pruppache and Klett 1997 eq.13-61
      libcloudphxx_const(si::dimensionless, f1, .78, 1)
      libcloudphxx_const(si::dimensionless, f2, .308, 1)

      //a, b coefficients from Morrison 2005 eq.4
      //tricky part is that the fall speed velocity parametrisation assumed here 
      //is based on droplet mass(in grams) Simmel et al. (2002) table 2
      //and not droplet diameter (in micro metres) Morrison 2005 eq.A4
      template<typename real_t>
      quantity<si::dimensionless, real_t> a_fall(
        const quantity<si::mass_density, real_t> &rhod,
        const quantity<si::mass_density, real_t> &rhod_rr,
        const quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> &rhod_nr
      ) {
        quantity<si::length, real_t> drop_r = r_drop_r(rhod_rr, rhod_nr);

        return rhod / rho_stp<real_t>()                        //to make it dimensionless         .... kilograms to grams
               * alpha_fall(drop_r) * std::pow(c_md<real_t>() * si::cubic_metres / si::kilograms * 1000, beta_fall(drop_r))
               * std::pow(real_t(1e-6), d_md<real_t>() * beta_fall(drop_r));
      }                         //^^^ metres to micro metres

      template<typename real_t>
      quantity<si::dimensionless, real_t> b_fall(
        const quantity<si::mass_density, real_t> &rhod_rr,
        const quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> &rhod_nr
      ) {
        quantity<si::length, real_t> drop_r = r_drop_r(rhod_rr, rhod_nr);

        return d_md<real_t>() * beta_fall(drop_r);
      }

      // relaxation time for condensation/evaporation term for rain drops
      // (see eq.4 in Morrison 2005)
      template<typename real_t>
      quantity<si::time, real_t> tau_relax_r(
        const quantity<si::temperature, real_t> &T, 
        const quantity<si::mass_density, real_t> &rhod,
        const quantity<si::mass_density, real_t> &rhod_rr,
        const quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> &rhod_nr
      ) {
        return real_t(1) / (
          (real_t(2) * pi<real_t>() * D_0<real_t>() * N0_r(rhod_nr, rhod_rr) * std::tgamma(real_t(2)))
	   * (
	     f1<real_t>() 
	     / (std::pow(lambda_r(rhod_nr, rhod_rr) * si::metres, real_t(2)) / si::square_metres)
	     +
	     f2<real_t>() 
	     * std::sqrt(a_fall(rhod, rhod_rr, rhod_nr) * rhod / visc(T) * si::square_metres / si::seconds)
	     * std::pow(Sc(visc(T), rhod, D_0<real_t>()), real_t(1./3)) * std::tgamma((b_fall(rhod_rr, rhod_nr) + real_t(5)) / real_t(2.))
	     * std::pow(lambda_r(rhod_nr, rhod_rr) * si::metres, -(b_fall(rhod_rr, rhod_nr) + 5) / real_t(2.)) * si::square_metres
	   )
        );
      }

      //drv_s/dT (derived from Clapeyron equation and pv = rv * rho_d * R_v * T)
      typedef divide_typeof_helper<si::dimensionless, si::temperature>::type one_over_temperature;
      template<typename real_t> 
      quantity<one_over_temperature, real_t> drv_s_dT(
        const quantity<si::temperature, real_t> &T,
        const quantity<si::dimensionless, real_t> &r_vs
      ) {
        return l_v(T) * r_vs / R_v<real_t>() / (T*T);
      }

      //condensation/evaporation rate
      template<typename real_t>
      quantity<si::frequency, real_t> cond_evap_rate(
        const quantity<si::temperature, real_t> T, 
        const quantity<si::pressure, real_t> p,
        const quantity<si::dimensionless, real_t> r_v,
        const quantity<si::time, real_t> tau_relax
      ) {
        const quantity<si::dimensionless, real_t> _r_vs = r_vs(T,p);
        return (r_v - _r_vs) / tau_relax / (1 + drv_s_dT(T, _r_vs) * l_v(T) / c_p(r_v));
      }                                                                     //TODO check ^ is it c_p or c_p(r)

    };
  };
};
