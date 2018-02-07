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
      const arrinfo_t<real_t> courant_x,
      const arrinfo_t<real_t> courant_y,
      const arrinfo_t<real_t> courant_z,
      const std::map<enum chem_species_t, const arrinfo_t<real_t> > ambient_chem
    )
    {
      if (init_called)
        throw std::runtime_error("init() may be called just once");
      init_called = true;

      // sanity checks
      if (th.is_null() || rv.is_null() || rhod.is_null())
        throw std::runtime_error("passing th, rv and rhod is mandatory");

      // --------  init cell characteristics  --------
      // initialising Eulerian-Lagrandian coupling
      if (!courant_x.is_null() || !courant_y.is_null() || !courant_z.is_null())
      {
        if (n_dims == 0)
          throw std::runtime_error("Courant numbers passed in 0D setup");

        if (n_dims == 1 && (courant_x.is_null() || !courant_y.is_null() || !courant_z.is_null()))
          throw std::runtime_error("Only X Courant number allowed in 1D setup");

        if (n_dims == 2 && (courant_x.is_null() || !courant_y.is_null() || courant_z.is_null()))
          throw std::runtime_error("Only X and Z Courant numbers allowed in 2D setup");

        if (n_dims == 3 && (courant_x.is_null() || courant_y.is_null() || courant_z.is_null()))
          throw std::runtime_error("All XYZ Courant number components required in 3D setup");
      }

      if (opts_init.chem_switch && ambient_chem.size() != chem_gas_n)
        throw std::runtime_error("chemistry was not switched off and ambient_chem is empty");

      if (!opts_init.chem_switch && ambient_chem.size() != 0)
        throw std::runtime_error("chemistry was switched off and ambient_chem is not empty");

      // TODO: in source match chemical composition
      if (opts_init.chem_switch && opts_init.src_switch)
        throw std::runtime_error("chemistry and aerosol source are not compatible");

      if(opts_init.dry_distros.size() > 1 && opts_init.chem_switch)
        throw std::runtime_error("chemistry and multiple kappa distributions are not compatible");

      // TODO: in source match kappas 
      if(opts_init.dry_distros.size() > 1 && opts_init.src_switch)
        throw std::runtime_error("aerosol source and multiple kappa distributions are not compatible");

//      if(opts_init.dry_distros.size() > 0 && opts_init.dry_sizes.size() > 0)
  //      throw std::runtime_error("Define dry_distros or dry_sizes, not both");

      if(opts_init.dry_distros.size() == 0 && opts_init.dry_sizes.size() == 0)
        throw std::runtime_error("Both dry_distros and dry_sizes are undefined");

      if(opts_init.dry_sizes.size() > 1)
        throw std::runtime_error("Multiple kappas dry_sizes is not supported yet");

  //    if(opts_init.dry_sizes.size() > 0 && opts_init.sd_conc > 0)
//        throw std::runtime_error("Dry_sizes init is not compatible with sd_conc (works only with const_multi)");

      if(opts_init.dry_sizes.size() > 0 && opts_init.sd_const_multi_dry_sizes == 0)
        throw std::runtime_error("Dry_sizes init used. Please specify multiplicity of dry_sizes SDs: sd_const_multi_dry_sizes.");

      if(opts_init.dry_sizes.size() == 0 && opts_init.sd_const_multi_dry_sizes != 0)
        throw std::runtime_error("sd_const_multi_dry_sizes is only compatible with dry_sizes init, which is not used in this run (dry_sizes undefined)");

      if(opts_init.sd_conc_large_tail && opts_init.sd_conc == 0)
        throw std::runtime_error("Sd_conc_large_tail make sense only with sd_conc init (i.e. sd_conc>0)");

      if(opts_init.sd_const_multi > 0 && opts_init.src_switch)
        throw std::runtime_error("aerosol source and constant multiplicity option are not compatible");

        if (n_dims > 0)
        {
          if (!(opts_init.x0 >= 0 && opts_init.x0 < m1(opts_init.nx) * opts_init.dx))
            throw std::runtime_error("!(x0 >= 0 & x0 < min(1,nx)*dz)");
          if (!(opts_init.y0 >= 0 && opts_init.y0 < m1(opts_init.ny) * opts_init.dy))
            throw std::runtime_error("!(y0 >= 0 & y0 < min(1,ny)*dy)");
          if (!(opts_init.z0 >= 0 && opts_init.z0 < m1(opts_init.nz) * opts_init.dz))
            throw std::runtime_error("!(z0 >= 0 & z0 < min(1,nz)*dz)");
          // check temporarily disabled since dewv_id is not passed anymore, TODO: fix it
//        if (!(opts_init.x1 > opts_init.x0 && opts_init.x1 <= m1(opts_init.nx) * opts_init.dx) && dev_id == -1) // only for single device runs, since on multi_CUDA x1 is not yet adjusted to local domain
//            throw std::runtime_error("!(x1 > x0 & x1 <= min(1,nx)*dx)");
          if (!(opts_init.y1 > opts_init.y0 && opts_init.y1 <= m1(opts_init.ny) * opts_init.dy))
            throw std::runtime_error("!(y1 > y0 & y1 <= min(1,ny)*dy)");
          if (!(opts_init.z1 > opts_init.z0 && opts_init.z1 <= m1(opts_init.nz) * opts_init.dz))
            throw std::runtime_error("!(z1 > z0 & z1 <= min(1,nz)*dz)");
        }

        if (opts_init.dt == 0) throw std::runtime_error("please specify opts_init.dt");
        if (opts_init.sd_conc * opts_init.sd_const_multi != 0) throw std::runtime_error("specify either opts_init.sd_conc or opts_init.sd_const_multi, not both");
        if (opts_init.sd_conc == 0 && opts_init.sd_const_multi == 0 && opts_init.sd_const_multi_dry_sizes == 0) throw std::runtime_error("please specify opts_init.sd_conc, opts_init.sd_const_multi or opts_init.sd_const_multi_dry_sizes");
        if (opts_init.coal_switch)
        {
          if(opts_init.terminal_velocity == vt_t::undefined) throw std::runtime_error("please specify opts_init.terminal_velocity or turn off opts_init.coal_switch");
          if(opts_init.kernel == kernel_t::undefined) throw std::runtime_error("please specify opts_init.kernel");
        }
        if (opts_init.sedi_switch)
          if(opts_init.terminal_velocity == vt_t::undefined) throw std::runtime_error("please specify opts_init.terminal_velocity or turn off opts_init.sedi_switch");
    }
  };
};

