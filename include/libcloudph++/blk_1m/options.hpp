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
    struct opts_t {
      bool 
        cevp = true, 
        revp = true, 
        conv = true, 
        clct = true, 
        sedi = true;
      real_t dt = 0;
    };
//</listing>
  }; 
};
