/** @file
  * @copyright University of Warsaw
  * @brief Definition of a structure holding options for Lagrangian microphysics
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libcloudph++/lgrngn/extincl.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
//<listing>
    template<typename real_t>
    struct opts_t 
    {
      // process toggling
      bool adve, sedi, cond, coal; 

      // RH limit for drop growth
      real_t RH_max;       

      // no. of substeps 
      int sstp_cond, sstp_coal;    

      // ctor with defaults (C++03 compliant) ...
//</listing>
      opts_t() : 
        adve(1), sedi(1), cond(1), coal(1), // all on
        sstp_cond(10), sstp_coal(10),
        RH_max(44) // :) (anything greater than 1.1 would be enough
      {}
    };
  }
};
