/** @file
  * @copyright University of Warsaw
  * @brief Definition of a structure holding options for Lagrangian microphysics
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <boost/ptr_container/ptr_unordered_map.hpp>
#include "../common/unary_function.hpp"

namespace libcloudphxx
{
  namespace lgrngn
  {
    template<typename real_t>
    struct opts_t
    {
      // processes
      bool 
        adve = true, 
        cond = true, 
        sedi = true, 
        coal = true, 
        chem = false;
// TODO: vent?, recycling

      // initial dry spectra
      typedef boost::ptr_unordered_map<real_t, common::unary_function<real_t> > dry_distros_t;

      //
      int nx, ny, nz; 
      real_t sd_conc_mean; 
      real_t dx, dy, dz; 
      dry_distros_t dry_distros;

      //real_t dt = 0;
    };
  }
};
