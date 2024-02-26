// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief initialisation routine for super droplets
  */

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_sanity_check(
      const arrinfo_t<real_t> th,
      const arrinfo_t<real_t> rv,
      const arrinfo_t<real_t> rhod,
      const arrinfo_t<real_t> p,
      const arrinfo_t<real_t> courant_x,
      const arrinfo_t<real_t> courant_y,
      const arrinfo_t<real_t> courant_z,
      const std::map<enum chem_species_t, const arrinfo_t<real_t> > ambient_chem
    )
    {
      if (init_called)
        throw std::runtime_error("libcloudph++: init() may be called just once");
      init_called = true;

      // sanity checks
      if (th.is_null() || rv.is_null() || rhod.is_null())
        throw std::runtime_error("libcloudph++: passing th, rv and rhod is mandatory");


      // --------  init cell characteristics  --------
      // initialising Eulerian-Lagrandian coupling
      if (!courant_x.is_null() || !courant_y.is_null() || !courant_z.is_null())
      {
        if (n_dims == 0)
          throw std::runtime_error("libcloudph++: Courant numbers passed in 0D setup");

        if (n_dims == 1 && (courant_x.is_null() || !courant_y.is_null() || !courant_z.is_null()))
          throw std::runtime_error("libcloudph++: Only X Courant number allowed in 1D setup");

        if (n_dims == 2 && (courant_x.is_null() || !courant_y.is_null() || courant_z.is_null()))
          throw std::runtime_error("libcloudph++: Only X and Z Courant numbers allowed in 2D setup");

        if (n_dims == 3 && (courant_x.is_null() || courant_y.is_null() || courant_z.is_null()))
          throw std::runtime_error("libcloudph++: All XYZ Courant number components required in 3D setup");
      }

      if (opts_init.chem_switch && ambient_chem.size() != chem_gas_n)
        throw std::runtime_error("libcloudph++: chemistry was not switched off and ambient_chem is empty");

      if (!opts_init.chem_switch && ambient_chem.size() != 0)
        throw std::runtime_error("libcloudph++: chemistry was switched off and ambient_chem is not empty");

      if (opts_init.chem_switch && opts_init.src_type!=src_t::off)
        throw std::runtime_error("libcloudph++: chemistry and aerosol source are not compatible");

      if (opts_init.src_type!=src_t::off && opts_init.src_dry_distros.empty() && opts_init.src_dry_sizes.empty())
        throw std::runtime_error("libcloudph++: CCN source enabled, but src_dry_distros and src_dry_sizes are empty");

      if (opts_init.src_type!=src_t::off && opts_init.src_dry_distros.size() > 1)
        throw std::runtime_error("libcloudph++: src_dry_distros can only have a single kappa value.");

      if (opts_init.src_type==src_t::matching && opts_init.dry_distros.size() > 1)
        throw std::runtime_error("libcloudph++: For 'matching' CCN source, the initial aerosol distribution can only have one kappa value (na kappa matching done).");

      if (opts_init.src_type!=src_t::off && n_dims<2)
        throw std::runtime_error("libcloudph++: CCN source works in 2D and 3D only.");

      if (opts_init.src_type==src_t::matching && !opts_init.src_dry_distros.empty() &&
          opts_init.src_dry_distros.begin()->first != opts_init.dry_distros.begin()->first) throw std::runtime_error("libcloudph++: For 'matching' CCN source, kappa of the source has to be the same as that of the initial profile (no kappa matching done)");

      if(opts_init.dry_distros.size() > 1 && opts_init.chem_switch)
        throw std::runtime_error("libcloudph++: chemistry and multiple kappa distributions are not compatible");

      if(opts_init.dry_distros.size() == 0 && opts_init.dry_sizes.size() == 0)
        throw std::runtime_error("libcloudph++: Both dry_distros and dry_sizes are undefined");

      if(opts_init.sd_conc_large_tail && opts_init.sd_conc == 0)
        throw std::runtime_error("libcloudph++: Sd_conc_large_tail make sense only with sd_conc init (i.e. sd_conc>0)");

      if(opts_init.sd_const_multi > 0 && opts_init.src_type!=src_t::off)
        throw std::runtime_error("libcloudph++: aerosol source and constant multiplicity option are not compatible"); // NOTE: why not?

      if (n_dims > 0)
      {
        if (!(opts_init.x0 >= 0 && opts_init.x0 < m1(opts_init.nx) * opts_init.dx))
          throw std::runtime_error("libcloudph++: !(x0 >= 0 & x0 < min(1,nx)*dz)");
        if (!(opts_init.y0 >= 0 && opts_init.y0 < m1(opts_init.ny) * opts_init.dy))
          throw std::runtime_error("libcloudph++: !(y0 >= 0 & y0 < min(1,ny)*dy)");
        if (!(opts_init.z0 >= 0 && opts_init.z0 < m1(opts_init.nz) * opts_init.dz))
          throw std::runtime_error("libcloudph++: !(z0 >= 0 & z0 < min(1,nz)*dz)");
        // check temporarily disabled since dewv_id is not passed anymore, TODO: fix it
//      if (!(opts_init.x1 > opts_init.x0 && opts_init.x1 <= m1(opts_init.nx) * opts_init.dx) && dev_id == -1) // only for single device runs, since on multi_CUDA x1 is not yet adjusted to local domain
//          throw std::runtime_error("libcloudph++: !(x1 > x0 & x1 <= min(1,nx)*dx)");
        if (!(opts_init.y1 > opts_init.y0 && opts_init.y1 <= m1(opts_init.ny) * opts_init.dy))
          throw std::runtime_error("libcloudph++: !(y1 > y0 & y1 <= min(1,ny)*dy)");
        if (!(opts_init.z1 > opts_init.z0 && opts_init.z1 <= m1(opts_init.nz) * opts_init.dz))
          throw std::runtime_error("libcloudph++: !(z1 > z0 & z1 <= min(1,nz)*dz)");
      }

      if (opts_init.dt == 0) throw std::runtime_error("libcloudph++: please specify opts_init.dt");
      if (opts_init.sd_conc * opts_init.sd_const_multi != 0) throw std::runtime_error("libcloudph++: specify either opts_init.sd_conc or opts_init.sd_const_multi, not both");
      if (opts_init.sd_conc == 0 && opts_init.sd_const_multi == 0 && opts_init.dry_sizes.size() == 0) throw std::runtime_error("libcloudph++: please specify opts_init.sd_conc, opts_init.sd_const_multi or opts_init.dry_sizes");
      if (opts_init.coal_switch)
      {
        if(opts_init.terminal_velocity == vt_t::undefined) throw std::runtime_error("libcloudph++: please specify opts_init.terminal_velocity or turn off opts_init.coal_switch");
        if(opts_init.kernel == kernel_t::undefined) throw std::runtime_error("libcloudph++: please specify opts_init.kernel");
      }
      if (opts_init.sedi_switch)
        if(opts_init.terminal_velocity == vt_t::undefined) throw std::runtime_error("libcloudph++: please specify opts_init.terminal_velocity or turn off opts_init.sedi_switch");
      if (opts_init.sedi_switch && opts_init.nz == 0)
        throw std::runtime_error("libcloudph++: opts_init.sedi_switch can be True only if n_dims > 1");
      if (opts_init.subs_switch && opts_init.nz == 0)
        throw std::runtime_error("libcloudph++: opts_init.subs_switch can be True only if n_dims > 1");
      if (opts_init.turb_adve_switch && opts_init.nz == 0)
        throw std::runtime_error("libcloudph++: opts_init.turb_adve_switch can be True only if n_dims > 1");
      if (opts_init.turb_cond_switch && opts_init.nz == 0)
        throw std::runtime_error("libcloudph++: opts_init.turb_cond_switch can be True only if n_dims > 1");
      if (opts_init.subs_switch && opts_init.nz != w_LS.size())
        throw std::runtime_error("libcloudph++: opts_init.subs_switch == True, but subsidence velocity profile size != nz");
      if ((opts_init.turb_adve_switch || opts_init.turb_cond_switch) && opts_init.nz != SGS_mix_len.size())
        throw std::runtime_error("libcloudph++: at least one of opts_init.turb_adve_switch, opts_init.turb_cond_switch is true, but SGS mixing length profile size != nz");
      if(opts_init.SGS_mix_len.size() > 0 && *std::min(opts_init.SGS_mix_len.begin(), opts_init.SGS_mix_len.end()) <= 0)
        throw std::runtime_error("libcloudph++: SGS_mix_len <= 0");
      if (!opts_init.aerosol_conc_factor.empty() && n_dims<2)
        throw std::runtime_error("libcloudph++: aerosol_conc_factor can only be used in 2D and 3D");
      if (!opts_init.aerosol_conc_factor.empty() && opts_init.nz != opts_init.aerosol_conc_factor.size())
        throw std::runtime_error("libcloudph++: aerosol_conc_factor size needs to be either 0 or nz");
      if (!opts_init.aerosol_conc_factor.empty() && opts_init.aerosol_independent_of_rhod==false)
      {
        std::cerr << "aerosol conc factor size: " << opts_init.aerosol_conc_factor.size() << std::endl;
        throw std::runtime_error("libcloudph++: aerosol_conc_factor can only be used if aerosol_independent_of_rhod==true");
        }
      #if defined(USE_MPI)
        if(opts_init.rlx_switch)
          std::cerr << "libcloudph++ WARNING: relaxation is not fully supported in MPI runs. Mean calculation and addition of SD will be done locally on each node." << std::endl;
        if(opts_init.chem_switch)
          throw std::runtime_error("libcloudph++: chemistry is not compatible with MPI");
      #endif
      if(n_dims < 2 && opts_init.rlx_switch)
        throw std::runtime_error("libcloudph++: CCN relaxation works only in 2D and 3D, set rlx_switch to false");
      if(opts_init.rlx_switch && opts_init.rlx_bins <= 0)
        throw std::runtime_error("libcloudph++: rlx_bins <= 0");
      if(opts_init.rlx_switch && opts_init.rlx_sd_per_bin <= 0)
        throw std::runtime_error("libcloudph++: rlx_sd_per_bin <= 0");
      if(opts_init.rlx_switch && opts_init.rlx_timescale <= 0)
        throw std::runtime_error("libcloudph++: rlx_timescale <= 0");        
      if(opts_init.rlx_switch && opts_init.chem_switch)
        throw std::runtime_error("libcloudph++: CCN relaxation does not work with chemistry");
    }
  };
};

