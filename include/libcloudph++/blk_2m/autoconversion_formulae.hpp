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
        real_t &rc,
        real_t &nc,
        const quantity<divide_typeof_helper<si::mass, si::volume>::type, real_t> &rhod
      ) {

        quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> N_c = rhod * nc / si::kilograms;
/*
std::cerr<<"autoconversion rate calc:"<<std::endl;
std::cerr<<"rc^2.47   = "<<  std::pow(rc, real_t(2.47)) <<std::endl; 
std::cerr<<"N_c^-1.79 = "<<  std::pow(N_c * real_t(1e-6) * si::cubic_metres, real_t(-1.79)) <<std::endl;
std::cerr<<"rhod      = "<<  rhod <<std::endl;
std::cerr<<"nc        = "<<  nc <<std::endl;
std::cerr<<"N_c       = "<<  N_c <<std::endl;
*/
        return real_t(1350) / si::seconds
          * std::pow(rc, real_t(2.47)) 
          * std::pow(N_c * real_t(1e-6) * si::cubic_metres, real_t(-1.79));
        //                        ^^^^
        //                        \__  m-3  -->  cm-3 
      }

    };
  };
};
