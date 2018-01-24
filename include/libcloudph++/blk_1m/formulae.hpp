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
      //libcloudphxx_const(si::frequency, k_autoconv_default, .001, si::hertz)

      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> autoconversion_rate(
	const quantity<si::dimensionless, real_t> &rc,
	const quantity<si::dimensionless, real_t> rc_thresh,
	//const quantity<si::frequency, real_t> k_autoconv = k_autoconv_default<real_t>()
	const quantity<si::frequency, real_t> k_autoconv
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
	return k_2<real_t>() * rc * std::pow(rr, real_t(.875));
      }

      // Kessler evaporation rate
      // eq. 5c in Grabowski & Smolarkiewicz 1996 (multiplied by rho!)
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> evaporation_rate( 
	quantity<si::dimensionless, real_t> rv,
	quantity<si::dimensionless, real_t> rvs,
	quantity<si::dimensionless, real_t> rr,
	quantity<si::mass_density, real_t> rhod,
	quantity<si::pressure, real_t> p
      )
      {
	return 
          (1 - rv / rvs) / rhod
	  * (
            real_t(1.6) 
            + real_t(124.9) * std::pow( 
              real_t(1e-3) * rhod * rr * si::cubic_metres / si::kilograms,
              real_t(.2046)
            ) 
          ) // ventilation factor TODO- move to ventil.hpp
	  * std::pow(
            real_t(1e-3) * rhod * rr * si::cubic_metres / si::kilograms, 
            real_t(.525)
          ) 
	  / (real_t(5.4e2) 
          + real_t(2.55e5) 
          * (real_t(1) / (p / si::pascals) / rvs)) 
	  / si::seconds * si::kilograms / si::cubic_metres;
      }

      // Kessler/Beard terminal velocity
      // eq. 5d in Grabowski & Smolarkiewicz 1996	
      libcloudphxx_const(si::velocity, vterm_A, 36.34, si::metre_per_second)

      using inverse_density = divide_typeof_helper<si::dimensionless,si::mass_density>::type;
      libcloudphxx_const(inverse_density, vterm_B, 1e-3, si::cubic_metres / si::kilograms)

      template <typename real_t>
      quantity<si::velocity, real_t> v_term(
	const quantity<si::dimensionless, real_t> &rr,
	const quantity<si::mass_density, real_t> &rhod,
	const quantity<si::mass_density, real_t> &rhod_0
      ) {
	return 
          vterm_A<real_t>() 
          * real_t(std::pow(
            (rhod * rr * vterm_B<real_t>()), 
            real_t(.1346)
          ) 
          * sqrt(rhod_0 / rhod)
        );
      }
    };
  };
};
