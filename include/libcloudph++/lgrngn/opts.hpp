/** @file
  * @copyright University of Warsaw
  * @brief Definition of a structure holding options for Lagrangian microphysics
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libcloudph++/lgrngn/extincl.hpp>
#include <libcloudph++/lgrngn/chem.hpp>

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
//</listing>

      // chem stuff
      bool chem;

      int sstp_chem; 

      std::vector<real_t> chem_gas;

      // ctor with defaults (C++03 compliant) ...
      opts_t() : 
        adve(true), sedi(true), cond(true), coal(true), chem(true), // all on 
        sstp_cond(10), sstp_coal(10), sstp_chem(10),
        RH_max(44), // :) (anything greater than 1.1 would be enough
        chem_gas(chem_gas_n)
      {
        for(int i=0; i<chem_gas_n; ++i)
        {
          chem_gas[i] = 0;
        }
      }
    };
  }
};
