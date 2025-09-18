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
      // Khairoutdinov and Kogan 2000 (eq. 29)
      // (note: a density-based formula in SI units is available in Wood 2005, table 1) 
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

      // NOTE: only working combinations are: th_dry == true && const_p == false; th_dry == false && const_p == true
      bool th_dry  = true, // if true, input and output theta are dry-air potential temperature; if false, they are 'standard' potential temperature
           const_p = false; // if true, pressure is equal to a supplied profile except for solving velocity (e.g. anelastic model); if false, pressure comes from the gas equation
    };
//</listing>
  }
 };
