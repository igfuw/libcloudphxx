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
    struct opts
    {
      bool 
        acti = true, 
        cond = true, 
        accr = true, 
        acnv = true, 
        turb = true,
        sedi = true;
      real_t dt = 0;
    };
  }
};
