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
    void particles_t<real_t, device>::impl::init_SD_with_distros_tail(const common::unary_function<double> &fun, const real_t log_rd_min_init)
    {
      dist_analysis_const_multi(fun);
 
      // overwrite calculated log_rd_min with external one (equal to log_rd_max of sd_conc init)
      log_rd_min=log_rd_min_init;
      if(log_rd_min >= log_rd_max)
        throw std::runtime_error(detail::formatter() << "Distribution analysis error: rd_min(" << exp(log_rd_min) << ") >= rd_max(" << exp(log_rd_max) << ")");
      
      // init number of SDs of this kappa in cells, TODO: due to rounding, we might end up with not exactly sd_conc SDs per cell...
      init_count_num_const_multi(fun, 1); // const_multi=1 for tail prtcls
  
      // update no of particles
      // TODO: move to a separate function
      n_part_old = n_part;
      n_part_to_init = thrust::reduce(count_num.begin(), count_num.end());
      n_part += n_part_to_init;
      hskpng_resize_npart(); 
  
      // init ijk vector, also n_part and resize n_part vectors
      init_ijk();
  
      // initialising dry radii (needs ijk)
      init_dry_const_multi(fun);
  
      // init multiplicities
      init_n_const_multi(1); 
    }
  };
};
