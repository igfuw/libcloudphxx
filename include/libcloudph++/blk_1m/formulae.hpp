/** @file
  * @copyright University of Warsaw
  * @brief single-moment bulk parameterisation formulae (Kessler)
  *   from @copybrief bib::Grabowski_and_Smolarkiewicz_1996
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libcloudph++/common/moist_air.hpp>

namespace libcloudphxx
{
  namespace blk_1m
  {
    namespace formulae
    {
      //Kessler autoconversion
      //eq. 5a in Grabowski & Smolarkiewicz 1996
      libcloudphxx_const(si::frequency, k_autoconv_default, .001, si::hertz)
      libcloudphxx_const(si::dimensionless, rc_thresh_default, .0005, 1) 

      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> autoconversion_rate(
	const quantity<si::dimensionless, real_t> &rc,
	const quantity<si::dimensionless, real_t> rc_thresh = rc_thresh_default<real_t>(),
	const quantity<si::frequency, real_t> k_autoconv = k_autoconv_default<real_t>()
      ) {
	return k_autoconv * std::max( real_t(0) * si::dimensionless(), rc - rc_thresh);
      }

      //Kessler collection
      //eq. 5b in Grabowski & Smolarkiewicz 1996
      libcloudphxx_const(si::frequency, k_2, 2.2, si::hertz)

      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> collection_rate(
	const quantity<si::dimensionless, real_t> &rc,
	const quantity<si::dimensionless, real_t> &rr
      ) {
	return k_2<real_t>() * rc * pow(rr, .875);
      }

      // Kessler evaporation rate
      // eq. 5c in Grabowski & Smolarkiewicz 1996 (multiplied by rho!)
      template <typename real_t>
      quantity<divide_typeof_helper<si::mass_density, si::time>::type, real_t> evaporation_rate( 
	quantity<si::dimensionless, real_t> rv,
	quantity<si::dimensionless, real_t> rvs,
	quantity<si::mass_density, real_t> rhod_rr,
	quantity<si::pressure, real_t> p
      )
      {
	return (1 - rv / rvs) 
	  * (1.6 + 124.9 * pow(1e-3 * (rhod_rr * si::cubic_metres / si::kilograms), .2046)) // ventilation factor TODO- move to ventil.hpp
	  * pow(1e-3 * (rhod_rr * si::cubic_metres / si::kilograms), .525) 
	  / (5.4e2 + 2.55e5 * (1. / (p / si::pascals) / rvs)) 
	  / si::seconds * si::kilograms / si::cubic_metres;
      }

      // Kessler/Beard terminal velocity
      // eq. 5d in Grabowski & Smolarkiewicz 1996	
      libcloudphxx_const(si::velocity, vterm_A, 36.34, si::metre_per_second)

      using inverse_density = divide_typeof_helper<si::dimensionless,si::mass_density>::type;
      libcloudphxx_const(inverse_density, vterm_B, 1e-3, si::cubic_metres / si::kilograms)

      template <typename real_t>
      quantity<si::velocity, real_t> v_term(
	const quantity<si::mass_density, real_t> &rho_r,
	const quantity<si::mass_density, real_t> &rho_d,
	const quantity<si::mass_density, real_t> &rho_d0
      ) {
	return vterm_A<real_t>() * real_t(pow(rho_r * vterm_B<real_t>(), .1346) * sqrt(rho_d0 / rho_d));
      }
    };
  };
};
