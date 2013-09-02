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

      //eq.2 Morrison and Grabowski 2007
      template<typename real_t>
      quantity<si::dimensionless, real_t> eta(
        quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> N
      ) {
        return .0005714 * 1e-6 * N * si::cubic_metres + .2714;
      }                 //^^^ convert N to 1/cm3

      //assumed mass-diametrer relationship coefficients (below A2 and A3 in Morrison 2005)
      //(mass of droplet = volume * density of water)
      libcloudphxx_const_derived(si::mass_density, c_md, pi<real_t>()/6 * rho_w<real_t>())
      libcloudphxx_const(si::dimensionless,        d_md, 3, 1)

      //for cloud droplets gamma size distribution is assumed 
      //(notation as in Morrison and Grabowski 2007 eq.1)

      //spectral index
      template<typename real_t>
      quantity<si::dimensionless, real_t> miu_c(
        quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> N
      ) {
        return 1 / pow<2>(eta(N)) - 1;
      }
      //slope
      template<typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::length>::type, real_t> lambda_c(
         quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> N,
         quantity<divide_typeof_helper<si::mass, si::volume>::type, real_t> rhod_rc
      ) {
        return pow(
                   c_md<real_t>() * N * tgamma(miu_c(N) + d_md<real_t>() + 1) / (rhod_rc * tgamma(miu_c(N) + 1)) * si::cubic_metres, 
                   real_t(1) / d_md<real_t>()
                  ) / si::metres;
      }
      //intercept  
      template<typename real_t>
      quantity<power_typeof_helper<si::length, static_rational<-4>>::type, real_t> N0_c(
         quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> N,
         quantity<divide_typeof_helper<si::mass, si::volume>::type, real_t> rhod_rc
      ) {
        return N * pow(lambda_c(N, rhod_rc) * si::metres, miu_c(N) + 1) / tgamma(miu_c(N) + 1) / si::metres;
      }

      //for rain drops Marshal-Palmer (exponential) size distribution is assumed
      //(meaning in fact that the spectral index from gamma distribtion miu is assumed to be 0)
      //slope
      template<typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::length>::type, real_t> lambda_r(
         quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> N,
         quantity<divide_typeof_helper<si::mass, si::volume>::type, real_t> rhod_rr
      ) {
        return pow(
                   c_md<real_t>() * N * tgamma(d_md<real_t>() +1) / rhod_rr * si::cubic_metres, 
                   real_t(1) / d_md<real_t>()
                  ) / si::metres;
      }
      //intercept  
      template<typename real_t>
      quantity<power_typeof_helper<si::length, static_rational<-4>>::type, real_t> N0_r(
         quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> N,
         quantity<divide_typeof_helper<si::mass, si::volume>::type, real_t> rhod_rr
      ) {
        return N * lambda_c(N, rhod_rr);
      }
      
      //cloud droplet radius based on rhod_rc, N (mean of gamma size distribution)
      template<typename real_t>
      quantity<si::length, real_t> r_drop_c(
        quantity<divide_typeof_helper<si::mass, si::volume>::type, real_t> rhod_rc,
        quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> N
      ) {
        if (rhod_rc > 0 * si::kilograms / si::cubic_metres && N > 0 / si::cubic_metres) 
          {return (miu_c(N) + 1) / lambda_c(N, rhod_rc) / real_t(2);}
        else {return 0 * si::metres;}
      } 

      //rain drop radius based on rhod_rr, N (mean of exponential size dist.)
      template<typename real_t>
      quantity<si::length, real_t> r_drop_r(
        quantity<divide_typeof_helper<si::mass, si::volume>::type, real_t> rhod_rr,
        quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> N
      ) {
        if (rhod_rr > 0 * si::kilograms / si::cubic_metres && N > 0 / si::cubic_metres)
          {return real_t(1) / lambda_r(N, rhod_rr) / real_t(2);}
        else {return 0 * si::metres;}
      } 

    };
  };
};
