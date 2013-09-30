/** @file
  * @copyright University of Warsaw
  * @brief Definition of a structure holding options for Lagrangian microphysics
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <boost/ptr_container/ptr_unordered_map.hpp>
#include <vector>

#include "../common/unary_function.hpp"
#include "../common/units.hpp"

namespace libcloudphxx
{
  namespace lgrngn
  {
//<listing>
    template<typename real_t>
    struct opts_t // C++03 compliant
    {
      // type definitions:
      typedef boost::ptr_unordered_map<
        // kappa
        real_t, 
        // n(ln(rd)) @ STP 
        common::unary_function<real_t> 
      > dry_distros_t;

      // member fields:
      bool 
        adve, sedi, cond, coal; // processes
      int 
        nx, ny, nz,             // grid
        sstp_cond, sstp_coal;   // substeps 
      real_t 
        dx, dy, dz, dt,         // steps
        sd_conc_mean, // super-droplet conc.       
        RH_max;       // condens. RH cutoff 
      dry_distros_t 
        dry_distros;  // initial dry aerosol 

      // methods ...
//</listing>
      // constructor with default values
      opts_t() : 
        adve(1), sedi(1), cond(1), coal(1),
        nx(0), ny(0), nz(0), 
        sstp_cond(10), sstp_coal(10),
        dx(1), dy(1), dz(1), 
        sd_conc_mean(64),
        RH_max(1.02) 
      {}
    };
  }
};
