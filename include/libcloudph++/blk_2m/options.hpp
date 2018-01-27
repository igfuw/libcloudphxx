/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief Definition of a structure holding options for double-moment bulk microphysics 
  */

#pragma once

#include <vector>

namespace libcloudphxx
{
  namespace blk_2m
  {
//<listing>
    template<typename real_t>
    struct opts_t
    {
      bool 
        acti = true, // activation
        cond = true, // condensation
        acnv = true, // autoconversion
        accr = true, // accretion
        sedi = true; // sedimentation

      // RH limit for activation
      real_t RH_max = 44; 
 
      // default parameters in the autoconversion parameterization
      real_t 
        acnv_A = real_t(1350),
        acnv_b = real_t(2.47),
        acnv_c = real_t(-1.79);
     
      // aerosol spectrum 
      struct lognormal_mode_t 
      { 
        real_t
          mean_rd,   // [m]
          sdev_rd,   // [1]
          N_stp,     // [m-3] @STP
          chem_b;    // [1]
      };
      std::vector<lognormal_mode_t> dry_distros;
    };
//</listing>
  }
 };
