// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */
// #include <limits>
#include <thrust/unique.h>
#include <thrust/binary_search.h>

namespace libcloudphxx
{
  namespace lgrngn
  {
    using std::get;

    // create new aerosol particles based on a size distribution
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::src_dry_distros_simple(const src_dry_distros_t<real_t> &sdd)
    {   
      // We assume that sdd size is 1
      // TODO: add a loop to allow sdd.size>1
      const auto p_sdd = sdd.begin()

      // add the source only once every number of steps
      if(src_stp_ctr % get<2>(p_sdd->second) != 0) return;

      const real_t sup_dt = get<2>(p_sdd->second) * opts_init.dt;

      // set number of SDs to init; use count_num as storage
      init_count_num_src(get<1>(p_sdd->second));

      // analyze distribution to get rd_min and max needed for bin sizes
      // TODO: this could be done once at the beginning of the simulation
      dist_analysis_sd_conc(
        *(get<0>(p_sdd->second)),
        get<1>(p_sdd->second),
        sup_dt
      ); 

      namespace arg = thrust::placeholders;

      // set no of particles to init
      n_part_old = n_part;
      n_part_to_init = thrust::reduce(count_num.begin(), count_num.end());
      n_part = n_part_old + n_part_to_init;
      hskpng_resize_npart();

      thrust_size_t n_part_bfr_src = n_part_old,
                    n_part_tot_in_src = n_part_to_init;

      // init ijk and rd3 of new particles
      init_ijk();
      init_dry_sd_conc(); 

      // init other peoperties of SDs that didnt have a match
      init_kappa(
        p_sdd->first.first
      ); 
      init_ice(
        p_sdd->first.second
      ); 

      if(opts_init.diag_incloud_time)
        init_incloud_time();

      init_n_sd_conc(
        *get<0>(p_sdd->second)
      ); 

      // init rw
      init_wet();
    
      // ijk -> i, j, k
      unravel_ijk(n_part_old);

      // init x, y, z, i, j, k
      init_xyz();

      // TODO: init chem
    }
  };  
};
