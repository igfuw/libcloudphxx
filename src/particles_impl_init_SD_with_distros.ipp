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
      for (typename opts_init_t<real_t>::dry_distros_t::const_iterator ddi = opts_init.dry_distros.begin(); ddi != opts_init.dry_distros.end(); ++ddi)
      {
        if(opts_init.sd_conc > 0)
          dist_analysis_sd_conc(
            ddi->second,
            opts_init.sd_conc
          );
/*
        else if(opts_init.sd_const_multi > 0)
          dist_analysis_const_multi(
            ddi->second
          );
*/
        tot_lnrd_rng += log_rd_max - log_rd_min;
      }

      // initialize SDs of each kappa-type
      for (typename opts_init_t<real_t>::dry_distros_t::const_iterator ddi = opts_init.dry_distros.begin(); ddi != opts_init.dry_distros.end(); ++ddi)
      {
        // another iteration needed to add SDs representing the large tail
        // they are added using the const_multi=1 logic
        const int n_iters = opts_init.sd_conc_large_tail ? 2 : 1; 

        real_t log_rd_max_sd_conc; // to store log_rd_max after sd_conc dist analysis (needed by tail)

        bool tail = 0; // are we initialising additional large tail prtcls
        for(int iter=0; iter < n_iters; ++iter, tail=!tail)
        { 
          // analyze the distribution, TODO: just did it
          if(opts_init.sd_conc > 0 && !tail)
          {
            dist_analysis_sd_conc(
              ddi->second,
              opts_init.sd_conc
            );
            log_rd_max_sd_conc = log_rd_max;
          }
          else if(opts_init.sd_const_multi > 0)
            dist_analysis_const_multi(
              ddi->second
            );
          else if(tail)
          {
            dist_analysis_const_multi(
              ddi->second
            );
            // tail particles will be initialised between log_rd_max of sd_conc and log_rd_max of const_multi
            log_rd_min = log_rd_max_sd_conc; 
          }
          if(log_rd_min >= log_rd_max)
            throw std::runtime_error("Distribution analysis error: rd_min >= rd_max");
          
          // init number of SDs of this kappa in cells, TODO: due to rounding, we might end up with not exactly sd_conc SDs per cell...
          if(opts_init.sd_conc > 0 && !tail)
          {
            // adjust the multiplicity init coefficient to smaller number of SDs representing this kappa-type
            real_t fraction = (log_rd_max - log_rd_min) / tot_lnrd_rng;
            multiplier *= opts_init.sd_conc / int(fraction * opts_init.sd_conc + 0.5);
            init_count_num_sd_conc(fraction);
          }
          else if(opts_init.sd_const_multi > 0)
            init_count_num_const_multi(ddi->second);
          else if(tail)
            init_count_num_const_multi(ddi->second, 1); // const_multi=1 for tail prtcls
  
          // update no of particles
          // TODO: move to a separate function
          n_part_old = n_part;
          n_part_to_init = thrust::reduce(count_num.begin(), count_num.end());
          n_part += n_part_to_init;
          hskpng_resize_npart(); 

          init_sstp();
  
          // init ijk vector, also n_part and resize n_part vectors
          init_ijk();
  
          // initialising dry radii (needs ijk)
          if(opts_init.sd_conc > 0 && !tail)
            init_dry_sd_conc();
          else if(opts_init.sd_const_multi > 0 || tail)
            init_dry_const_multi(ddi->second);

          // init kappa
          init_kappa(ddi->first);
  
          // init multiplicities
          if(opts_init.sd_conc > 0 && !tail)
            init_n_sd_conc(ddi->second); // TODO: document that n_of_lnrd_stp is expected!
          else if(opts_init.sd_const_multi > 0)
            init_n_const_multi(); 
          else if(tail)
            init_n_const_multi(1); 
  
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
      }
    }
  };
};
