/** @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date July 2013
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include <libcloudph++/common/units.hpp> // TODO: do detail?
#include <libcloudph++/common/macros.hpp> // TODO: do detail?
#include <libcloudph++/common/moist_air.hpp> // TODO: do detail?

namespace libcloudphxx
{
  namespace common
  {
    namespace kelvin
    {
      //done to avoid BOOST_COMMA in preprocesor macro below
      typedef divide_typeof_helper<si::force, si::length>::type force_over_length;
      // water surface tension (Petters and Kreidenweis 2007 //TODO - check)
      libcloudphxx_const(force_over_length, sg_surf, 0.072, si::newton/si::metres)

      // Kelvin curvature parameter (see eq. 7 in Kvorostyanov and Curry 2006)
      template <typename real_t>
      quantity<si::length, real_t> A(
        quantity<si::temperature, real_t> T
      ) {
        return real_t(2.) * sg_surf<real_t>() / moist_air::R_v<real_t>() / T / moist_air::rho_w<real_t>();
      }

    }
  }
};

