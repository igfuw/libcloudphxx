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
      
      //droplet radius based on rhod_rc and concentration
      template<typename real_t>
      quantity<si::length, real_t> r_drop(
        quantity<divide_typeof_helper<si::mass, si::volume>::type, real_t> rhod_rc,
        quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> N
      ) {
        return root<3>(3. * rhod_rc / (4. * pi<real_t>() * N * rho_w<real_t>()));
      } 

    };
  };
};
