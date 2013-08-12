/** @file
  * @copyright University of Warsaw
  * @brief double-moment bulk condensation/evaporation parameterisation formulae
  * @date August 2013
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once
#include <libcloudph++/common/moist_air.hpp>
#include <libcloudph++/blk_2m/common_formulae.hpp>

namespace libcloudphxx
{
  namespace blk_2m
  {
    namespace formulae
    { 
      //terminal fall velocity based on the data from Gunn and Kinzer (1949) and Beard (1976) 
      //modified by Simmel et al. (2002) -> see table 2 there

      libcloudphxx_const(si::dimensionless, a1, 4.5795 * 1e5, 1)
      libcloudphxx_const(si::dimensionless, a2, 4.962  * 1e3, 1)
      libcloudphxx_const(si::dimensionless, a3, 1.732  * 1e3, 1)
      libcloudphxx_const(si::dimensionless, a4, 9.17   * 1e2, 1)

      libcloudphxx_const(si::dimensionless, b1, 2./3, 1)
      libcloudphxx_const(si::dimensionless, b2, 1./3, 1)
      libcloudphxx_const(si::dimensionless, b3, 1./6, 1)
      libcloudphxx_const(si::dimensionless, b4, 0   , 1)

      libcloudphxx_const(si::length, d1, 134.43  * 1e-6, si::metres)
      libcloudphxx_const(si::length, d2, 1511.64 * 1e-6, si::metres)
      libcloudphxx_const(si::length, d3, 3477.84 * 1e-6, si::metres)

      template <typename real_t>
      quantity<si::velocity, real_t> v_term(
        const quantity<si::mass_density, real_t> &rhod_rr,
        const quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> &rhod_nr
      ) {                                                           
          if (real_t(2) * r_drop(rhod_rr, rhod_nr) < d1<real_t>())  //.... mass of drop in grams
            return a1<real_t>() * pow(rhod_rr/rhod_nr/si::kilograms * 1000, b1<real_t>()) * real_t(1e-2) * si::metres/si::seconds;
          //                                          ^^^^^^^^^^^^^ to make it dimensionless       ^^^^ velocity in metres/seconds

          else if (real_t(2) * r_drop(rhod_rr, rhod_nr) < d2<real_t>())
            return a2<real_t>() * pow(rhod_rr/rhod_nr/si::kilograms * 1000, b2<real_t>()) * real_t(1e-2) * si::metres/si::seconds;

          else if (real_t(2) * r_drop(rhod_rr, rhod_nr) < d3<real_t>())
            return a3<real_t>() * pow(rhod_rr/rhod_nr/si::kilograms * 1000, b3<real_t>()) * real_t(1e-2) * si::metres/si::seconds;

          else 
            return a4<real_t>() * pow(rhod_rr/rhod_nr/si::kilograms * 1000, b4<real_t>()) * real_t(1e-2) * si::metres/si::seconds;
      }
    };
  };
};
