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
    { //autoconversion rate
      
      //as in Khairoutdinov and Kogan 2000
      //but taken from Wood 2005 (table 1) - si units!

      libcloudphxx_const(si::dimensionless, A_auto, 7.42 * 1e13, 1)

      //in activation term all activated droples are assumed to have the radius of 1 um
      libcloudphxx_const(si::length, drizzle_radius, 25 * 1e-6, si::metres);

      template<typename real_t>
      quantity<si::frequency, real_t> autoconv_rate(
        quantity<si::mass_density, real_t> rhod,
        quantity<si::mass_density, real_t> rhod_rc,
        quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> N_c
      ) {
        return A_auto<real_t>() / si::seconds
               * pow(rhod_rc * si::cubic_metres / si::kilograms, 2.47) 
               * pow(N_c * si::cubic_metres, -1.79) 
               * pow(rhod * si::cubic_metres / si::kilograms, -1.47);
      }

    };
  };
};
