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

      // water surface tension (as used in Petters and Kreidenweis 2007)
      libcloudphxx_const(force_over_length, sg_surf, 0.072, si::newton/si::metres)

      // Kelvin curvature parameter (see eq. 7 in Kvorostyanov and Curry 2006)
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::length, real_t> A(
        quantity<si::temperature, real_t> T
      ) {
        using namespace moist_air;
        return real_t(2) * sg_surf<real_t>() / R_v<real_t>() / T / rho_w<real_t>();
      }

      // Kelvin term in Koehler equation (see eq. 1 in Petters and Kreidenweis 2007)
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::dimensionless, real_t> klvntrm(
	quantity<si::length, real_t> r,
	quantity<si::temperature, real_t> T
      )   
      {   
	return exp(A<real_t>(T) / r);
      }   
    }
  }
};

