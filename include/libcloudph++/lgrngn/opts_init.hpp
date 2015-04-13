/** @file
  * @copyright University of Warsaw
  * @brief Definition of a structure holding options for Lagrangian microphysics
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libcloudph++/lgrngn/extincl.hpp>
#include <libcloudph++/lgrngn/kernel.hpp>
#include <libcloudph++/lgrngn/terminal_velocity.hpp>
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

      // no. of substeps 
      int sstp_cond, sstp_coal; 

      // Lagrangian domain extents
      real_t x0, y0, z0, x1, y1, z1;

      // mean no. of super-droplets per cell
      real_t sd_conc_mean; 

      // coalescence Kernel type
      kernel_t kernel;

      // terminal velocity formula
      vt_t terminal_velocity;
//</listing>
 
      // coalescence kernel parameters
      std::vector<real_t> kernel_parameters;

      // chem
      bool chem_switch;  // if false no chemical reactions throughout the whole simulation (no memory allocation)
      int sstp_chem;
      real_t chem_rho;

      // RH threshold for calculating equilibrium condition at t=0
      real_t RH_max;

      // ctor with defaults (C++03 compliant) ...
      opts_init_t() : 
        nx(0), ny(0), nz(0), // the defaults are OK for a parcel set-up 
        dx(1), dy(1), dz(1), // (but are only used to compute n_part -
        x0(0), y0(0), z0(0), //  dv is computed from rhod assuming 
        x1(1), y1(1), z1(1), //  that the parcel contains 1kg of dry air)
        sd_conc_mean(0), 
        dt(0),   
        sstp_cond(1), sstp_coal(1), sstp_chem(1),         
        chem_switch(false),  // chemical reactions turned off by default
        RH_max(.95), // value seggested in Lebo and Seinfeld 2011
        chem_rho(0), // dry particle density  //TODO add checking if the user gave a different value (np w init)  (was 1.8e-3)
        terminal_velocity(undefined)
      {}
    };
  }
};
