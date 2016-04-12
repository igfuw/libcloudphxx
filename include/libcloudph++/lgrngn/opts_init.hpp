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
  
      // timestep interval at which source will be applied
      int supstp_src;

      // Lagrangian domain extents
      real_t x0, y0, z0, x1, y1, z1;

      // no. of super-droplets per cell
      unsigned long long sd_conc; 

      // max no. of super-droplets in the system
      // should be enough to store particles from sources
      unsigned long long n_sd_max; 

      // source distro per unit time
      dry_distros_t src_dry_distros;

      // number of SDs created per cell per source iteration
      unsigned long long src_sd_conc;

      // height up to which aerosol will be created
      // will be rounded to cell number - cells are supposed to be uniform
      real_t src_z1;

      // coalescence Kernel type
      kernel_t::kernel_t kernel;

      // terminal velocity formula
      vt_t::vt_t terminal_velocity;
//</listing>
 
      // coalescence kernel parameters
      std::vector<real_t> kernel_parameters;

      // chem
      bool chem_switch,  // if false no chemical reactions throughout the whole simulation (no memory allocation)
           coal_switch,  // if false no coalescence throughout the whole simulation
           sedi_switch,  // if false no sedimentation throughout the whole simulation
           src_switch;  // if false no source throughout the whole simulation

      int sstp_chem;
      real_t chem_rho;

      // RH threshold for calculating equilibrium condition at t=0
      real_t RH_max;

      // rng seed
      int rng_seed;

      // no of GPUs per MPI node to use, 0 for all available
      int dev_count; 

      // ctor with defaults (C++03 compliant) ...
      opts_init_t() : 
        nx(0), ny(0), nz(0),
        dx(1), dy(1), dz(1),
        x0(0), y0(0), z0(0),
        x1(1), y1(1), z1(1),
        sd_conc(0), 
        dt(0),   
        sstp_cond(1), sstp_coal(1), sstp_chem(1),         
        supstp_src(1),
        chem_switch(false),  // chemical reactions turned off by default
        sedi_switch(true),  // sedimentation turned on by default
        coal_switch(true),  // coalescence turned on by default
        src_switch(false),  // source turned off by default
        RH_max(.95), // value seggested in Lebo and Seinfeld 2011
        chem_rho(0), // dry particle density  //TODO add checking if the user gave a different value (np w init)  (was 1.8e-3)
        rng_seed(44),
        terminal_velocity(vt_t::undefined),
        kernel(kernel_t::undefined),
        dev_count(0),
        n_sd_max(0),
        src_sd_conc(0),
        src_z1(0)
      {}

      // dtor (just to silence -Winline warnings)
      ~opts_init_t() {}
    };
  }
};
