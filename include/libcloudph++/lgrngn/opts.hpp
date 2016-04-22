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
      bool adve, sedi, cond, coal, src, rcyc;

      // RH limit for drop growth
      real_t RH_max;       
//</listing>

      // process toggling for chemistry
      bool chem_dsl, chem_dsc, chem_rct;
      // toggling between open and closed chemical system
      bool chem_sys_cls;

      // ctor with defaults (C++03 compliant) ...
      opts_t() : 
        adve(true), sedi(true), cond(true), coal(true), src(false), rcyc(false),
        chem_dsl(false), chem_dsc(false), chem_rct(false), chem_sys_cls(false),
        RH_max(44) // :) (anything greater than 1.1 would be enough
      {
      }
    };
  }
};
