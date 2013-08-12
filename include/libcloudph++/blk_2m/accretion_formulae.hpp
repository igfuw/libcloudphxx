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
    { //accretion rate
      
      //as in Khairoutdinov and Kogan 2000
      //but taken from Wood 2005 (table 1) - (all hail to si units!)

      libcloudphxx_const(si::dimensionless, A_acc, 67, 1)

      template<typename real_t>
      quantity<divide_typeof_helper<si::mass_density, si::time>::type, real_t> accretion_rate(
        quantity<si::mass_density, real_t> rhod,
        quantity<si::mass_density, real_t> rhod_rc,
        quantity<si::mass_density, real_t> rhod_rr
      ) {
        return A_acc<real_t>() * si::kilograms / si::cubic_metres / si::seconds
               * pow(rhod_rc * si::cubic_metres / si::kilograms * rhod_rr * si::cubic_metres / si::kilograms, 1.15) 
               * pow(rhod * si::cubic_metres / si::kilograms, -1.3);
      }

    };
  };
};
