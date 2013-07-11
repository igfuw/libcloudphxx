/** @file
  * @copyright University of Warsaw
  * @brief Definition of a structure holding options for Lagrangian microphysics
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

namespace libcloudphxx
{
  namespace lgrngn
  {
    template<typename real_t>
    struct opts
    {
      bool 
        adve = true, 
        cond = true, 
        sedi = true, 
        coal = true, 
        chem = false;
// TODO: vent?
      real_t dt = 0;
    };
  }
};
