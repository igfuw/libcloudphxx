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
      auto p_sdd = sdd.cbegin();

      if(std::get<1>(p_sdd->first) < 0 || std::get<1>(p_sdd->first) > 1)
        throw std::runtime_error("libcloudph++: soluble_fraction in opts.src_dry_distros must be in [0, 1]");

      const int sd_conc = get<2>(p_sdd->first);
      const int sd_const_multi = get<3>(p_sdd->first);
      const int supstp = get<4>(p_sdd->first);
      if(sd_conc > 0 && sd_const_multi > 0)
        throw std::runtime_error("libcloudph++: specify either sd_conc or sd_const_multi for each source dry distribution, not both");

      // add the source only once every number of steps
      assert(supstp > 0);
      if(src_stp_ctr % supstp != 0) return;

      const real_t sup_dt = supstp * opts_init.dt;

      if(sd_conc > 0)
      {
        init_count_num_src(sd_conc);
        init_dist_analysis_sd_conc(*p_sdd->second, sd_conc, sup_dt);
      }
      else if(sd_const_multi > 0)
      {
        init_dist_analysis_const_multi(*p_sdd->second);
        init_count_num_src_const_multi(*p_sdd->second, sd_const_multi, sup_dt);
      }
      else
        throw std::runtime_error("libcloudph++: specify either sd_conc or sd_const_multi for each source dry distribution");

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
      if(sd_conc > 0)
      {
        init_dry_sd_conc();
        init_n_sd_conc(*p_sdd->second);
      }
      else
      {
        init_dry_const_multi(*p_sdd->second);
        init_n_const_multi(sd_const_multi);
      }

      // init other properties of SDs
      init_kappa(
        std::get<0>(p_sdd->first),
        std::get<1>(p_sdd->first)
      );
      if (opts_init.ice_switch)
      {
        init_insol(std::get<1>(p_sdd->first));
        init_a_c_rho_ice();
        if (! opts_init.time_dep_ice_nucl)
        {
          init_T_freeze();
        }
      }

      if(opts_init.diag_incloud_time)
        init_incloud_time();

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
