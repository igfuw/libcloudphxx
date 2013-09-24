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
//<listing>
    template<typename real_t>
    struct opts_t    // uses C++11
    {
      bool 
        acti = true, // activation
        cond = true, // condensation
        acnv = true, // autoconversion
        accr = true, // accretion
        sedi = true; // sedimentation
      
      // aerosol 
      real_t
        mean_rd, // [m]
        sdev_rd, // [1]
        N_stp,   // [m-3] @STP
        chem_b;  // [1]
    };
//</listing>
  }
 };
