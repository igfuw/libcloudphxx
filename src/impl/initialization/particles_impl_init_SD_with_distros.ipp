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
    // init SD parameters from a dry size distribution
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_SD_with_distros()
    {
      // initialize SDs of each kappa-type
      for (auto ddi = opts_init.dry_distros.cbegin(); ddi != opts_init.dry_distros.cend(); ++ddi)
      {
        const auto &distro = ddi->second;
        const auto sd_conc = std::get<2>(ddi->first);
        const auto sd_const_multi = std::get<3>(ddi->first);
        if(sd_conc > 0)
        {
          init_SD_with_distros_sd_conc(*distro, sd_conc);
          init_SD_with_distros_finalize(ddi->first);
          
          if(opts_init.sd_conc_large_tail)
          {
            init_SD_with_distros_tail(*distro, log_rd_max);
            init_SD_with_distros_finalize(ddi->first);
          }
        }
        if(sd_const_multi > 0)
        {
          init_SD_with_distros_const_multi(*distro, sd_const_multi);
          init_SD_with_distros_finalize(ddi->first);
        }
      }
    }

    // final inits common for tail/sd_conc/const_multi
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_SD_with_distros_finalize(const std::tuple<real_t, real_t, unsigned long long, unsigned long long> &kpa_sol_frac, const bool unravel_ijk_switch)
    {
      // dry_distros defines total dry radius; insoluble part determined by soluble_fraction

      // init kappa
      init_kappa(std::get<0>(kpa_sol_frac), std::get<1>(kpa_sol_frac));

      if (opts_init.ice_switch)
      {
        init_insol(std::get<1>(kpa_sol_frac));

        init_a_c_rho_ice();
        if (! opts_init.time_dep_ice_nucl)
        {
          init_T_freeze();
        }
      }
      
      // initialising wet radii
      init_wet();

      // memory allocation for chemical reactions, done after init.grid to have npart defined
      if(opts_init.chem_switch){
        init_chem();
      }

      // initialising mass of chemical compounds in droplets (needs to be done after dry radius)
      if(opts_init.chem_switch){
        init_chem_aq();
      }
      
      // init for substepping for chem reactions
      if(opts_init.chem_switch){
       init_percell_sstp_chem();
      }

      // calculate initail volume (helper for Henry in chem)
      if (opts_init.chem_switch){
        chem_vol_ante();
      }
      
      // ijk -> i, j, k
      if(unravel_ijk_switch)
        unravel_ijk(n_part_old);

      // initialising particle positions
      init_xyz();

      if(opts_init.diag_incloud_time)
        init_incloud_time();
    }
  };
};
