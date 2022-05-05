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
    // create new aerosol particles based on a size distribution
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::src_dry_distros_simple(const real_t &dt)
    {   
      // set number of SDs to init; use count_num as storage
      init_count_num_src(opts_init.src_sd_conc);

      // analyze distribution to get rd_min and max needed for bin sizes
      // TODO: this could be done once at the beginning of the simulation
      dist_analysis_sd_conc(
        *(opts_init.src_dry_distros.begin()->second),
        opts_init.src_sd_conc,
        dt
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
        opts_init.src_dry_distros.begin()->first
      ); 

      if(opts_init.diag_incloud_time)
        init_incloud_time();

      init_n_sd_conc(
        *(opts_init.src_dry_distros.begin()->second)
      ); // TODO: document that n_of_lnrd_stp is expected!
    
      // ijk -> i, j, k
      unravel_ijk(n_part_old);

      // init x, y, z
      init_xyz();

      hskpng_ijk_ref(n_part_old);

      // init rw
      init_wet();

      // TODO: init chem
    }
  };  
};
