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

      // eq.2 Morrison and Grabowski 2007
      template<typename real_t>
      inline quantity<si::dimensionless, real_t> eta(
        const quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> &n
      ) {
        return real_t(.0005714 * 1e-6) * n * si::cubic_metres + real_t(.2714);
      }                        //^^^ convert N to 1/cm3 

      // assumed mass-diametrer relationship coefficients (below A2 and A3 in Morrison 2005)
      // (mass of droplet = volume * density of water)
      libcloudphxx_const_derived(si::mass_density, c_md, pi<real_t>()/6 * rho_w<real_t>())
      libcloudphxx_const(si::dimensionless,        d_md, 3, 1)

      // for cloud droplets gamma size distribution is assumed 
      // (notation as in Morrison and Grabowski 2007 eq.1)

      // spectral index
      template<typename real_t>
      inline quantity<si::dimensionless, real_t> miu_c(
        const quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> &n
      ) {
        auto tmp = real_t(1) / pow<2>(eta(n)) - real_t(1);
        assert(finite(tmp) && "spectral index n is finite failed");
        return tmp;
      }

      // slope
      template<typename real_t>
      inline quantity<divide_typeof_helper<si::dimensionless, si::length>::type, real_t> lambda_c(
         const quantity<divide_typeof_helper<si::dimensionless, si::mass>::type, real_t> &n,
         const quantity<si::dimensionless, real_t> &rc,
        const quantity<si::mass_density, real_t> &rhod
      ) {
        auto tmp = pow(
	  c_md<real_t>() * n * std::tgamma(miu_c(n*rhod) + d_md<real_t>() + real_t(1)) 
          / 
	  (rc * std::tgamma(miu_c(n*rhod) + real_t(1))) * si::cubic_metres
	  , 
	  real_t(1) / d_md<real_t>()
        ) / si::metres;
        assert(finite(tmp * si::metres) && "slope lambda_c is finite failed");
        return tmp;
      }

      // intercept  
      template<typename real_t>
      inline quantity<divide_typeof_helper<
        power_typeof_helper<si::length, static_rational<-3>>::type,
        si::mass
      >, real_t> N0_c(
         const quantity<divide_typeof_helper<si::dimensionless, si::mass>::type, real_t> &nc,
         const quantity<si::dimensionless, real_t> &rc,
         const quantity<si::mass_density, real_t> &rhod
      ) {
        auto tmp = nc 
          * pow(lambda_c(nc, rc, rhod) * si::metres, miu_c(nc * rhod) + real_t(1)) 
          / std::tgamma(miu_c(nc * rhod) + real_t(1)) 
          / si::metres;
        assert(finite(tmp * si::metres * si::cubic_metres) && "intercept N0_c is finite failed");
        return tmp;
      }

      // for rain drops Marshal-Palmer (exponential) size distribution is assumed
      // (meaning in fact that the spectral index from gamma distribtion miu is assumed to be 0)
      // slope
      template<typename real_t>
      inline quantity<divide_typeof_helper<si::dimensionless, si::length>::type, real_t> lambda_r(
         const quantity<divide_typeof_helper<si::dimensionless, si::mass>::type, real_t> &nr,
         const quantity<si::dimensionless, real_t> &rr
      ) {
        auto tmp = pow(
	  c_md<real_t>() * nr * std::tgamma(d_md<real_t>() + real_t(1)) / rr * si::cubic_metres
          , 
	  real_t(1) / d_md<real_t>()
	) / si::metres;
        assert(finite(tmp * si::metres) && "slope lambda_r is finite failed");
        return tmp;
      }

      // intercept  
      template<typename real_t>
      inline quantity<divide_typeof_helper<
        power_typeof_helper<si::length, static_rational<-1>>::type,
        si::mass
      >::type, real_t> N0_r(
         const quantity<divide_typeof_helper<si::dimensionless, si::mass>::type, real_t> &n,
         const quantity<si::dimensionless, real_t> &rr
      ) {
        auto tmp = n * lambda_r(n, rr);
        assert(finite(tmp * si::kilograms * si::metres) && "intercept N0_r is finite failed");
        return tmp;
      }
      
      // cloud droplet radius based on rc, N (mean of gamma size distribution)
      template<typename real_t>
      inline quantity<si::length, real_t> r_drop_c(
        real_t &rc,
        real_t &nc,
        const quantity<si::mass_density, real_t> &rhod
      ) {
        quantity<divide_typeof_helper<si::dimensionless, si::mass>::type, real_t> n_c = nc / si::kilograms;
        quantity<si::dimensionless, real_t> r_c = rc * si::dimensionless();
        if (rc > 0 && nc > 0) 
          return (miu_c(nc / si::kilograms * rhod) + real_t(1)) / lambda_c(n_c, r_c, rhod) / real_t(2);
        else 
          return 0 * si::metres;
      } 

      // rain drop radius based on rr, N (mean of exponential size dist.)
      template<typename real_t>
      inline quantity<si::length, real_t> r_drop_r(
        const quantity<si::dimensionless, real_t> &rr,
        const quantity<divide_typeof_helper<si::dimensionless, si::mass>::type, real_t> &nr
      ) {
        if (rr > 0 && nr * si::kilograms > 0)
          return real_t(1) / lambda_r(nr, rr) / real_t(2);
        else 
          return 0 * si::metres;
      } 
    };
  };
};
