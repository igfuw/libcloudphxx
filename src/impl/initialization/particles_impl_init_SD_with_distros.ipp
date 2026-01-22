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
      // calc sum of ln(rd) ranges of all distributions
      real_t tot_lnrd_rng = 0.;
      if(opts_init.sd_conc > 0)
        for (auto ddi = opts_init.dry_distros.cbegin(); ddi != opts_init.dry_distros.cend(); ++ddi)
        {
            init_dist_analysis_sd_conc(
              *(ddi->second),
              opts_init.sd_conc
            );
          tot_lnrd_rng += log_rd_max - log_rd_min;
        }

      // initialize SDs of each kappa-type
      for (auto ddi = opts_init.dry_distros.cbegin(); ddi != opts_init.dry_distros.cend(); ++ddi)
      {
        if(opts_init.sd_conc > 0)
        {
          init_SD_with_distros_sd_conc(*(ddi->second), tot_lnrd_rng);
          init_SD_with_distros_finalize(ddi->first);
          
          if(opts_init.sd_conc_large_tail)
          {
            init_SD_with_distros_tail(*(ddi->second), log_rd_max);
            init_SD_with_distros_finalize(ddi->first);
          }
        }
        if(opts_init.sd_const_multi > 0)
        {
          init_SD_with_distros_const_multi(*(ddi->second));
          init_SD_with_distros_finalize(ddi->first);
        }
      }
    }

    // final inits common for tail/sd_conc/const_multi
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_SD_with_distros_finalize(const kappa_rd_insol_t<real_t> &kpa_rd_insol, const bool unravel_ijk_switch)
    {
      // init kappa
      init_kappa(kpa_rd_insol.kappa);

      if (opts_init.ice_switch)
      {
        init_insol_dry_sizes(kpa_rd_insol.rd_insol);
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
