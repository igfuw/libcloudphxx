/** @file
  * @copyright University of Warsaw
  * @brief Definition of a structure holding options for Lagrangian microphysics
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include "extincl.hpp"
#include "kernel.hpp"
#include "terminal_velocity.hpp"
#include "advection_scheme.hpp"
#include "RH_formula.hpp"
#include "../common/chem.hpp"

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
      // defined with a distribution
      // uses shared_ptr to make opts_init copyable
      typedef std::unordered_map<
        real_t,                // kappa
        std::shared_ptr<unary_function<real_t>> // n(ln(rd)) @ STP; alternatively it's n(ln(rd)) independent of rhod if aerosol_independent_of_rhod=true
      > dry_distros_t;
      dry_distros_t dry_distros;

      // defined with a size-number pair
      typedef std::map<
        real_t,                  // kappa
        std::map<real_t,         // radius [m]
          std::pair<real_t, int> // STP_concentration [1/m^3], SD multiplicity
        > 
      > dry_sizes_t;
      dry_sizes_t dry_sizes;

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
 
      // should more SDs be added to better represent large tail of the distribution
      bool sd_conc_large_tail;

      // should aerosol concentration init be independent of rhod (assumed to be in cm^{-3} and not at STP)
      bool aerosol_independent_of_rhod;

      // or, alternatively to sd_conc_mean, multiplicity of all SDs = const
      unsigned long long  sd_const_multi;

      // max no. of super-droplets in the system
      // should be enough to store particles from sources
      unsigned long long n_sd_max; 

      // source distro per unit time
      dry_distros_t src_dry_distros;

      // number of SDs created from src_dry_distros per cell per source iteration
      unsigned long long src_sd_conc;

      // dry sizes of droplets added from the source, STP_concentration created per unit time instead of the STP_concentration 
      dry_sizes_t src_dry_sizes;

      // box in which aerosol from source will be created
      // will be rounded to cell number - cells are supposed to be uniform
      real_t src_x0, src_y0, src_z0, src_x1, src_y1, src_z1;

      // coalescence Kernel type
      kernel_t::kernel_t kernel;

      // terminal velocity formula
      vt_t::vt_t terminal_velocity;

      // super-droplet advection scheme
      as_t::as_t adve_scheme;

      // RH formula
      RH_formula_t::RH_formula_t RH_formula;
//</listing>
 
      // coalescence kernel parameters
      std::vector<real_t> kernel_parameters;

      // chem
      bool chem_switch,  // if false no chemical reactions throughout the whole simulation (no memory allocation)
           coal_switch,  // if false no coalescence throughout the whole simulation
           sedi_switch,  // if false no sedimentation throughout the whole simulation
           subs_switch,  // if false no subsidence throughout the whole simulation
           src_switch,   // if false no source throughout the whole simulation
           exact_sstp_cond, // if true, use per-particle sstp_cond logic, if false, use per-cell
           turb_adve_switch,   // if true, turbulent motion of SDs is modeled
           turb_cond_switch,   // if true, turbulent condensation of SDs is modeled
           turb_coal_switch;   // if true, turbulent coalescence kernels can be used

      int sstp_chem;
      real_t chem_rho;

      // do we want to track the time SDs spend inside clouds
      bool diag_incloud_time;

      // RH threshold for calculating equilibrium condition at t=0
      real_t RH_max;

      // rng seed
      int rng_seed;

      // no of GPUs per MPI node to use, 0 for all available
      int dev_count; 

      // GPU number to use, only used in CUDA backend (and not in multi_CUDA)
      int dev_id;

      // subsidence rate profile, positive downwards [m/s]
      std::vector<real_t> w_LS;

      // SGS mixing length profile [m]
      std::vector<real_t> SGS_mix_len;

      real_t rd_min; // minimal dry radius of droplets (works only for init from spectrum)

      bool no_ccn_at_init; // if true, no ccn / SD are put at the start of the simulation

      bool open_side_walls,       // if true, side walls are "open", i.e. SD are removed at contact. Periodic otherwise.
           periodic_topbot_walls; // if true, top and bot walls are periodic. Open otherwise


      // ctor with defaults (C++03 compliant) ...
      opts_init_t() : 
        nx(0), ny(0), nz(0),
        dx(1), dy(1), dz(1),
        x0(0), y0(0), z0(0),
        x1(1), y1(1), z1(1),
        sd_conc(0), 
        sd_conc_large_tail(false), 
        aerosol_independent_of_rhod(false), 
        sd_const_multi(0),
        dt(0),   
        sstp_cond(1), sstp_coal(1), sstp_chem(1),         
        supstp_src(1),
        chem_switch(false),  // chemical reactions turned off by default
        sedi_switch(true),  // sedimentation turned on by default
        subs_switch(false),  // subsidence turned off by default
        coal_switch(true),  // coalescence turned on by default
        src_switch(false),  // source turned off by default
        exact_sstp_cond(false),
        turb_cond_switch(false),
        turb_adve_switch(false),
        turb_coal_switch(false),
        RH_max(.95), // value seggested in Lebo and Seinfeld 2011
        chem_rho(0), // dry particle density  //TODO add checking if the user gave a different value (np w init)  (was 1.8e-3)
        rng_seed(44),
        terminal_velocity(vt_t::undefined),
        kernel(kernel_t::undefined),
        adve_scheme(as_t::implicit),
        RH_formula(RH_formula_t::pv_cc),
        dev_count(0),
        dev_id(-1),
        n_sd_max(0),
        src_sd_conc(0),
        src_x0(0),
        src_x1(0),
        src_y0(0),
        src_y1(0),
        src_z0(0),
        src_z1(0),
        rd_min(0.),
        diag_incloud_time(false),
        no_ccn_at_init(false),
        open_side_walls(false),
        periodic_topbot_walls(false)
      {}

      // dtor (just to silence -Winline warnings)
      ~opts_init_t() {}
    };
  }
};
