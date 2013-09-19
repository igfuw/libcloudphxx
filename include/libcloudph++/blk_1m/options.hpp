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
    struct opts_t { // uses C++11
      bool 
        cevp = true, 
        revp = true, 
        conv = true, 
        clct = true, 
        sedi = true;
// TODO: autoconv threshold
    };
//</listing>
  }; 
};
