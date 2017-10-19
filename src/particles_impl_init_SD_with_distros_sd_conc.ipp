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
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_SD_with_distros_sd_conc(const std::shared_ptr<common::unary_function<real_t>> &fun, const real_t &tot_lnrd_rng)
    {
      // analyze the distribution, TODO: just did it in init_SD_with_distros
      dist_analysis_sd_conc(
          fun,
          opts_init.sd_conc
        );
      if(log_rd_min >= log_rd_max)
        throw std::runtime_error(detail::formatter() << "Distribution analysis error: rd_min(" << exp(log_rd_min) << ") >= rd_max(" << exp(log_rd_max) << ")");
      
      // init number of SDs of this kappa in cells, TODO: due to rounding, we might end up with not exactly sd_conc SDs per cell...
        // adjust the multiplicity init coefficient to smaller number of SDs representing this kappa-type
      real_t fraction = (log_rd_max - log_rd_min) / tot_lnrd_rng;
      multiplier *= opts_init.sd_conc / int(fraction * opts_init.sd_conc + 0.5);
      init_count_num_sd_conc(fraction);
  
      // update no of particles
      // TODO: move to a separate function
      n_part_old = n_part;
      n_part_to_init = thrust::reduce(count_num.begin(), count_num.end());
      n_part += n_part_to_init;
      hskpng_resize_npart(); 
  
      // init ijk vector, also n_part and resize n_part vectors
      init_ijk();
  
      // initialising dry radii (needs ijk)
      init_dry_sd_conc();

      // init multiplicities
      init_n_sd_conc(fun); // TODO: document that n_of_lnrd_stp is expected!
    }
  };
};
