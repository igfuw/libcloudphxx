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
        for (typename opts_init_t<real_t>::dry_distros_t::const_iterator ddi = opts_init.dry_distros.begin(); ddi != opts_init.dry_distros.end(); ++ddi)
        {
            dist_analysis_sd_conc(
              ddi->second,
              opts_init.sd_conc
            );
          tot_lnrd_rng += log_rd_max - log_rd_min;
        }

      // initialize SDs of each kappa-type
      for (typename opts_init_t<real_t>::dry_distros_t::const_iterator ddi = opts_init.dry_distros.begin(); ddi != opts_init.dry_distros.end(); ++ddi)
      {
        if(opts_init.sd_conc > 0)
        {
          init_SD_with_distros_sd_conc(ddi->second, tot_lnrd_rng);
          init_SD_with_distros_finalize(ddi->first);
          
          if(opts_init.sd_conc_large_tail)
          {
            init_SD_with_distros_tail(ddi->second, log_rd_max);
            init_SD_with_distros_finalize(ddi->first);
          }
        }
        if(opts_init.sd_const_multi > 0)
        {
          init_SD_with_distros_const_multi(ddi->second);
          init_SD_with_distros_finalize(ddi->first);
        }
      }
    }

    // final inits common for tail/sd_conc/const_multi
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_SD_with_distros_finalize(const real_t &kappa)
    {
      // init kappa
      init_kappa(kappa);

      // save initial th/rv/rhod for substepping
      init_sstp();
      
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
       init_sstp_chem();
      }

      // calculate initail volume (helper for Henry in chem)
      if (opts_init.chem_switch){
        chem_vol_ante();
      }
      
      // initialising particle positions
      init_xyz();
    }
  };
};
