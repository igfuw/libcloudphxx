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
    void particles_t<real_t, device>::impl::src_dry_sizes(const real_t &dt)
    {
      using dry_sizes_t = typename opts_init_t<real_t>::dry_sizes_t;
      using kappa_t  = typename dry_sizes_t::key_type;
      using size_number_t = typename dry_sizes_t::mapped_type;
      //using conc_multi_t = typename size_number_t::mapped_type;


      // loop over kappas
      for (typename dry_sizes_t::const_iterator dsi = opts_init.src_dry_sizes.begin(); dsi != opts_init.src_dry_sizes.end(); ++dsi)
      {
        const kappa_t &kappa(dsi->first);
        const size_number_t &size_number_map(dsi->second);

        // loop over the "size : {concentration per second, multiplicity}" pairs for this kappa
        for (typename size_number_t::const_iterator sni = size_number_map.begin(); sni != size_number_map.end(); ++sni)
        {
          // init number of SDs of this kappa in cells
          init_count_num_src(sni->second.second);
  
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

          // init kappa
          init_kappa(kappa);
  
          // init multiplicities
          init_n_dry_sizes(sni->second.first*dt, sni->second.second); 
  
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
