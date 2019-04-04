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
    { // autoconversion rate

      // as in Khairoutdinov and Kogan 2000 (eq. 29)
      // (note: a density-based formula in SI units is available in Wood 2005, table 1)

      // in autoconversion term all drizzle drops are assumed to have the radius of 25 um
      libcloudphxx_const(si::length, drizzle_radius, 25 * 1e-6, si::metres);

      template<typename real_t>
      inline quantity<si::frequency, real_t> autoconv_rate(
        const real_t &rc,
        const real_t &nc,
        const quantity<divide_typeof_helper<si::mass, si::volume>::type, real_t> &rhod,
        const quantity<si::dimensionless, real_t> acnv_A,
        const quantity<si::dimensionless, real_t> acnv_b,
        const quantity<si::dimensionless, real_t> acnv_c
      ) {

        quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> N_c = rhod * nc / si::kilograms;

        return acnv_A / si::seconds
          * std::pow(rc, acnv_b)
          * std::pow(N_c * real_t(1e-6) * si::cubic_metres, acnv_c);
        //                        ^^^^
        //                        \__  m-3  -->  cm-3
      }

    };
  };
};
