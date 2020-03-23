/** @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date May 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once

#include "units.hpp"
#include "macros.hpp"

namespace libcloudphxx
{
  namespace common
  {
    namespace lognormal
    {
#if !defined(__NVCC__)
      using std::pow;
      using std::exp;
      using std::sqrt;
#endif
      // lognormal distribution as a function of ln(r) (Seinfeld & Pandis 1997 eq 7.33)
      template <typename real_t>
      quantity<power_typeof_helper<si::length,static_rational<-3>>::type, real_t> n_e(
	const quantity<si::length, real_t> &mean_r,
	const quantity<si::dimensionless, real_t> &stdev, 
	const quantity<power_typeof_helper<si::length,static_rational<-3>>::type, real_t> &n_tot, 
	const quantity<si::dimensionless, real_t> &lnr
      ) {
	return n_tot 
	  * real_t(exp(-pow((lnr - log(mean_r/si::metres)), 2) / real_t(2) / pow(log(stdev),2)))
	  / real_t(log(stdev))
	  / real_t(sqrt(2*pi<real_t>()))
	;
      }

      // lognormal distribution as a function of r (Seinfeld & Pandis 1997 eq 7.34)
      template <typename real_t>
      quantity<power_typeof_helper<si::length,static_rational<-4>>::type, real_t> n(
	const quantity<si::length, real_t> &mean_r,
	const quantity<si::dimensionless, real_t> &stdev, 
	const quantity<power_typeof_helper<si::length,static_rational<-3>>::type, real_t> &n_tot, 
	const quantity<si::length, real_t> &r
      ) {
	return n_tot / r
	  * real_t(exp(-pow((log(r/mean_r)), 2) / real_t(2) / pow(log(stdev),2)))
	  / real_t(log(stdev))
	  / real_t(sqrt(2*pi<real_t>()))
	;
      }
    };
  };
};
