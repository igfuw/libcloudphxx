/** @file
  * @copyright University of Warsaw
  * @brief Definition of a structure holding options for single-moment bulk microphysics
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

namespace libcloudphxx
{
  namespace blk_1m 
  { 
//<listing>
    template<typename real_t>
    struct opts_t {  // uses C++11
      bool 
        cond = true, // condensation
        cevp = true, // evaporation of cloud
        revp = true, // evaporation of rain 
        conv = true, // autoconversion
        accr = true, // accretion
        sedi = true; // sedimentation
      real_t 
        r_c0 = 5e-4; // autoconv. threshold
    };
//</listing>
  }; 
};
