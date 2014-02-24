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
    { // accretion rate
      
      // as in Khairoutdinov and Kogan 2000
      // (assumes that all the collected droplets have mean radius)
      template<typename real_t>
      inline quantity<divide_typeof_helper<si::frequency, si::mass>::type, real_t> collision_sink_rate(
        const quantity<si::frequency, real_t> &drr,
        const quantity<si::length, real_t> &r
      ) {
        return drr / (real_t(4./3.) * pi<real_t>() * (r*r*r) * common::moist_air::rho_w<real_t>());
      }

    };
  };
};
