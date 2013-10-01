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
    using common::unary_function;

    template<typename real_t>
    struct opts_init_t 
    {
      // initial dry sizes of aerosol
      typedef boost::ptr_unordered_map<
        real_t,                // kappa
        unary_function<real_t> // n(ln(rd)) @ STP 
      > dry_distros_t;
      dry_distros_t dry_distros;

      // Eulerian component parameters
      int nx, ny, nz;
      real_t dx, dy, dz, dt;

      // Lagrangian component parameters
      real_t
        sd_conc_mean; // super-droplets per cell

      // ctor with defaults (C++03 compliant) ...
      opts_init_t() : 
        nx(0), ny(0), nz(0), // parcel setup
        dx(1), dy(1), dz(1), // parcel setup
        sd_conc_mean(64)
      {}
    };

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
