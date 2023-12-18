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
    void particles_t<real_t, device>::impl::init_SD_with_sizes()
    {
//      using dry_sizes_t = typename opts_init_t<real_t>::dry_sizes_t;
  //    using size_number_t = typename dry_sizes_t::mapped_type;
      //using conc_multi_t = typename size_number_t::mapped_type;


      // loop over (kappa, ice) pairs
      for (auto dsi = opts_init.dry_sizes.cbegin(); dsi != opts_init.dry_sizes.cend(); ++dsi)
      {
        const real_t &kappa(dsi->first.first);
        const real_t &ice(dsi->first.second);
        const auto &size_number_map(dsi->second);

        // loop over the "size : {concentration, count}" pairs for this (kappa, ice) pair
        for (auto sni = size_number_map.cbegin(); sni != size_number_map.cend(); ++sni)
        {
          // init number of SDs of this kappa in cells
          init_count_num_dry_sizes(sni->second);
  
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
          init_ice(ice);
  
          // init multiplicities
          init_n_dry_sizes(sni->second.first, sni->second.second); 
  
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
        }
      }
    }
  };
};
