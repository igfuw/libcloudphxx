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
      // accretion rate as in Khairoutdinov and Kogan 2000
      // but taken from Wood 2005 (table 1) - (all hail to si units!)
      template<typename real_t>
      inline quantity<si::frequency, real_t> accretion_rate(
        const real_t &rc,
        const quantity<si::dimensionless, real_t> &rr
      ) {
        return real_t(67) / si::seconds * std::pow(rc * rr, real_t(1.15));
      }
    };
  };
};
