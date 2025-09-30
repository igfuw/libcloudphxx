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
    // initialize SD parameters with dry_radius-concentration pairs (i.e. dry_sizes)
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::src_dry_sizes(const src_dry_sizes_t<real_t> &sds)
    {
//      using dry_sizes_t = typename opts_t<real_t>::dry_sizes_t;
//      using real_t  = typename dry_sizes_t::key_type;
//      using size_number_t = typename dry_sizes_t::mapped_type;
      //using conc_multi_t = typename size_number_t::mapped_type;


      // loop over (kappa, ice) pairs
     // for (typename dry_sizes_t::const_iterator dsi = opts.src_dry_sizes.begin(); dsi != opts.src_dry_sizes.end(); ++dsi)
      for (auto dsi = sds.cbegin(); dsi != sds.cend(); ++dsi)
      {
        const real_t &kappa(dsi->first.first);
        const real_t &ice(dsi->first.second);
        const auto &size_number_map(dsi->second);

        // loop over the "size : {concentration per second, multiplicity, supstp}" for this (kappa, ice) pair
        for (auto sni = size_number_map.cbegin(); sni != size_number_map.cend(); ++sni)
        {
          // add the source only once every number of steps
          assert(get<2>(sni->second) > 0);
          if(src_stp_ctr % get<2>(sni->second) != 0) continue;

          const real_t sup_dt = get<2>(sni->second) * opts_init.dt;

          // init number of SDs of this kappa in cells
          init_count_num_src(get<1>(sni->second));
  
          // update no of particles
          // TODO: move to a separate function
          n_part_old = n_part;
          n_part_to_init = thrust::reduce(count_num.begin(), count_num.end());
          n_part += n_part_to_init;
          hskpng_resize_npart(); 

          // init ijk vector using count_num, also n_part and resize n_part vectors
          init_ijk();
  
          // initialising dry radii (needs ijk)
          init_dry_dry_sizes(sni->first);

          // init kappa and ice
          init_kappa(kappa);
          init_insol_dry_sizes(ice);
  
          // init multiplicities
          init_n_dry_sizes(get<0>(sni->second)*sup_dt, get<1>(sni->second)); 
  
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
      
          // ijk -> i, j, k
          unravel_ijk(n_part_old);
  
          // initialising particle positions
          init_xyz();

          if(opts_init.diag_incloud_time)
            init_incloud_time();
        }
      }
    }
  };
};
