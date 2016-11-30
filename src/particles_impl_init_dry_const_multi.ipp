// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief initialisation routine for super droplets
  */

#include <iostream>
#include <algorithm>

#include "detail/thrust.hpp"
#include "detail/functors_host.hpp"

#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/extrema.h>
#include <thrust/binary_search.h>

#include <libcloudph++/common/earth.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      // calculate cumulative distribution function
      template<typename real_t, typename vec_t>
      void calc_CDF(const common::unary_function<real_t> &fun, const real_t &min, const real_t &max, const real_t &bin_size, vec_t &vec)
      {
        const int n = (max - min) / bin_size + 1; //no of points at which cdf will be calculated
        vec.resize(n);

        namespace arg = thrust::placeholders;
        // fill vec with fun values at n points
        thrust::transform(
          thrust::make_transform_iterator(thrust::make_counting_iterator(0), min + bin_size * arg::_1),
          thrust::make_transform_iterator(thrust::make_counting_iterator(n), min + bin_size * arg::_1),
          vec.begin(), eval_and_mul<real_t>(fun, 1));

        // calculate CDF
        thrust::inclusive_scan(vec.begin(), vec.end(), vec.begin());
       
        // normalize CDF     
        thrust::transform(vec.begin(), vec.end(), vec.begin(), arg::_1 / vec.back());
      }
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_dry_const_multi(
      const common::unary_function<real_t> *n_of_lnrd_stp 
    )
    {
      // calculate cumulative distribution function
      thrust::host_vector<real_t> cdf;

      detail::calc_CDF(*n_of_lnrd_stp, log_rd_min, log_rd_max, config.bin_precision, cdf);

      // tossing random numbers [0,1] for dry radii
      rand_u01(n_part);

      // rd3 temporarily means logarithm of radius!
      thrust_device::vector<real_t> &lnrd(rd3);

      thrust::host_vector<real_t> host_u01(n_part); 
      thrust::copy(u01.begin(), u01.end(), host_u01.begin());
      thrust::host_vector<real_t> host_lnrd(n_part); 
      
      namespace arg = thrust::placeholders;
      // sample ln(rd) from the distribution with the inverse transform sampling method
      thrust::upper_bound(cdf.begin(), cdf.end(), host_u01.begin(), host_u01.end(), host_lnrd.begin());
      thrust::copy(host_lnrd.begin(), host_lnrd.end(), lnrd.begin());
      thrust::transform(lnrd.begin(), lnrd.end(), lnrd.begin(), log_rd_min + arg::_1 * config.bin_precision);

      // converting rd back from logarithms to rd3
      thrust::transform(
        lnrd.begin(),
        lnrd.end(),
        rd3.begin(),
        detail::exp3x<real_t>()
      );
    }
  };
};
