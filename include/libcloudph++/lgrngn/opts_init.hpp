/** @file
  * @copyright University of Warsaw
  * @brief Definition of a structure holding options for Lagrangian microphysics
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libcloudph++/lgrngn/extincl.hpp>
#include <libcloudph++/lgrngn/kernel.hpp>
#include <libcloudph++/lgrngn/chem.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    using common::unary_function;

//<listing>
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

      // Lagrangian domain extents
      real_t x0, y0, z0, x1, y1, z1;

      // mean no. of super-droplets per cell
      real_t sd_conc_mean; 

      // coalescence Kernel type
      kernel_t kernel;
//</listing>
      // chem
      real_t chem_rho;

      // ctor with defaults (C++03 compliant) ...
      opts_init_t() : 
        nx(0), ny(0), nz(0), // parcel setup, 1m3
        dx(1), dy(1), dz(1), // parcel setup, 1m3
        x0(0), y0(0), z0(0), // parcel setup, 1m3
        x1(1), y1(1), z1(1), // parcel setup, 1m3
        sd_conc_mean(64), 
        kernel(geometric),
        dt(1e-3),
        chem_rho(1.8e3) // dry particle density
      {}
    };
  }
};
