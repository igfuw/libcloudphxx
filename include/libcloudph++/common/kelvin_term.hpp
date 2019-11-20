/** @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date July 2013
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "units.hpp" // TODO: do detail?
#include "macros.hpp" // TODO: do detail?
#include "moist_air.hpp" // TODO: do detail?

namespace libcloudphxx
{
  namespace common
  {
    namespace kelvin
    {
      //done to avoid BOOST_COMMA in preprocesor macro below
      typedef divide_typeof_helper<si::force, si::length>::type force_over_length;

      //water - air surface tension, Eotvos rule, TODO: move somewhere else
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<force_over_length, real_t> sg_surf(
        quantity<si::temperature, real_t> T
      ) {
        return real_t(0.07275) * (real_t(1.) - real_t(0.002) * (T / si::kelvins - real_t(291.))) * si::newtons / si::metres;
      }

      // Kelvin curvature parameter (see eq. 7 in Kvorostyanov and Curry 2006)
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::length, real_t> A(
        quantity<si::temperature, real_t> T
      ) {
        using namespace moist_air;
        return real_t(2) * sg_surf<real_t>(T) / R_v<real_t>() / T / rho_w<real_t>();
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

