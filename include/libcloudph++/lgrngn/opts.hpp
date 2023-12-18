/** @file
  * @copyright University of Warsaw
  * @brief Definition of a structure holding options for Lagrangian microphysics
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include "extincl.hpp"
#include "distro_t.hpp"
#include "../common/chem.hpp"

namespace libcloudphxx
{
  namespace lgrngn
  {
//<listing>
    template<typename real_t>
    struct opts_t 
    {
      // process toggling
      bool adve, sedi, subs, cond, coal, src, rlx, rcyc, turb_adve, turb_cond, turb_coal;

      // RH limit for drop growth
      real_t RH_max;       
//</listing>

      // process toggling for chemistry
      bool chem_dsl, chem_dsc, chem_rct;

      // overriding dt from opts_init
      real_t dt;

      // aerosol source distro per unit time
      dry_distros_t<real_t> src_dry_distros;

      // dry sizes of droplets added from the source, STP_concentration created per unit time instead of the STP_concentration
      dry_sizes_t<real_t> src_dry_sizes;

      // ctor with defaults (C++03 compliant) ...
      opts_t() : 
        adve(true), sedi(true), subs(false), cond(true), coal(true), src(false), rlx(false), rcyc(false),
        chem_dsl(false), chem_dsc(false), chem_rct(false),
        turb_adve(false), turb_cond(false), turb_coal(false),
        RH_max(44), // :) (anything greater than 1.1 would be enough
        dt(-1) // negative means that we do not override dt in this step
      {
      }
    };
  }
};
