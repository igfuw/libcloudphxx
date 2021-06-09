// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief initialisation routine for super droplets
  */

#include <iostream>
#include <algorithm>

#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/extrema.h>

#include <libcloudph++/common/earth.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template<typename real_t>
      struct calc_lnrd
      {
        const real_t ln_rd_min,
                     ln_rd_max;
     
        calc_lnrd(const real_t &ln_rd_min, const real_t &ln_rd_max): ln_rd_min(ln_rd_min), ln_rd_max(ln_rd_max) {}
  
        template<typename Tuple>
        BOOST_GPU_ENABLED
        real_t operator()(Tuple tup)
        {
          return ln_rd_min + real_t(thrust::get<2>(tup) - thrust::get<3>(tup) + thrust::get<0>(tup)) * (ln_rd_max - ln_rd_min)  / real_t(thrust::get<1>(tup));
        }
      };
    };

    // init
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_dry_sd_conc()
    {
      // tossing random numbers
      rand_u01(n_part_to_init);

      // rd3 temporarily means logarithm of radius!
      thrust_device::vector<real_t> &lnrd(rd3);
      
      thrust_device::vector<thrust_size_t> &ptr(tmp_device_size_cell);
      thrust::exclusive_scan(count_num.begin(), count_num.end(), ptr.begin()); // number of SDs to init in cells up to (i-1), TODO: same is done in init_xyz, store it?
      
      // shifting from [0,1] to [log(rd_min),log(rd_max)] and storing into rd3
      // each log(radius) randomized only on a small subrange to make the distributions more uniform
      // particles are sorted by cell number (see init_ijk), uniform distribution in each cell
      // lnrd is not sorted
      thrust::transform(
        thrust::make_zip_iterator(thrust::make_tuple(
          u01.begin(),                                                              // random number
          thrust::make_permutation_iterator(count_num.begin(), ijk.begin() + n_part_old), // number of SDs in the cell
          zero,                                        // sequence to iterate over distribution
          thrust::make_permutation_iterator(ptr.begin(), ijk.begin() + n_part_old)        // number of SDs in cells up to this one
        )),
        thrust::make_zip_iterator(thrust::make_tuple(
          u01.begin(),                                                              // random number
          thrust::make_permutation_iterator(count_num.begin(), ijk.begin() + n_part_old), // number of SDs in the cell
          zero,                                        // sequence to iterate over distribution
          thrust::make_permutation_iterator(ptr.begin(), ijk.begin() + n_part_old)        // number of SDs in cells up to this one
        )) + n_part_to_init,
        lnrd.begin() + n_part_old, 
        detail::calc_lnrd<real_t>(log_rd_min, log_rd_max)
      );
      
      // converting rd back from logarithms to rd3
      thrust::transform(
        lnrd.begin() + n_part_old,
        lnrd.end(),
        rd3.begin() + n_part_old,
        detail::exp3x<real_t>()
      );
    }
  };
};
