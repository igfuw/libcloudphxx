/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief Definition of a structure holding options for double-moment bulk microphysics 
  */

#pragma once

namespace libcloudphxx
{
  namespace blk_2m
  {
    template<typename real_t>
    struct opts_t
    {
      bool 
        acti = true, 
        cond = true, 
        acnv = true, 
        accr = true, 
        sedi = true;
      
      real_t 
        dt = 0;

      //assumed aerosol size distribution (for activation)
      quantity<si::length, real_t> mean_rd;
      quantity<si::dimensionless, real_t> sdev_rd;
      quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t> N_tot;
      //assumed aerosol chemical composition (also for activation)
      quantity<si::dimensionless, real_t> chem_b;

    };
  }
 };
