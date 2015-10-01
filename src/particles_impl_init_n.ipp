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

#include <libcloudph++/common/earth.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t>
      struct eval_and_multiply
      {   
        const common::unary_function<real_t> &fun;
        const real_t &mul;

        // ctor
        eval_and_multiply(
          const common::unary_function<real_t> &fun, 
          const real_t &mul
        )
          : fun(fun), mul(mul)
        {}

        real_t operator()(real_t x)  // x is rd3
        {
          real_t lnrd = log(x) / 3.;
          return mul * fun(lnrd); 
        }
      };  
    };

    // init
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_n(
      const real_t kappa,
      const common::unary_function<real_t> *n_of_lnrd_stp // TODO: kappa-spectrum map
    )
    {
      // filling kappas
      thrust::fill(kpa.begin() + n_part_old, kpa.end(), kappa);

      // temporary space on the host 
      thrust::host_vector<real_t> tmp_real(n_part_to_init);
      thrust::host_vector<thrust_size_t> tmp_ijk(n_part_to_init);
      thrust::host_vector<real_t> &tmp_rhod(tmp_host_real_cell);

      thrust::copy(
	rhod.begin(), rhod.end(), // from
	tmp_rhod.begin()          // to
      );
      printf("rhod w init_n\n");
      debug::print(tmp_rhod);
      thrust::copy(
	ijk.begin()+n_part_old, ijk.end(), // from
	tmp_ijk.begin()         // to
      );
      printf("ijk w init_n\n");
      debug::print(tmp_ijk);

      // filling n with multiplicities
      // (performing it on a local copy as n_of_lnrd_stp may lack __device__ qualifier)
      // device -> host (not needed for omp or cpp ... but happens just once)
      thrust::copy(rd3.begin()+n_part_old, rd3.end(), tmp_real.begin());
      printf("rd3 w init_n\n");
      debug::print(tmp_real);
      
      // evaluating n_of_lnrd_stp
      thrust::transform(
        tmp_real.begin(), tmp_real.end(), // input 
        tmp_real.begin(),                 // output
        detail::eval_and_multiply<real_t>(*n_of_lnrd_stp, multiplier)
      );
      printf("n (real) po eval and multi w init_n, multiplier %lf\n", double(multiplier));
      debug::print(tmp_real);

      // correcting STP -> actual ambient conditions
      {
        namespace arg = thrust::placeholders;
        using common::earth::rho_stp;

        thrust::transform(
          tmp_real.begin(), tmp_real.end(),            // input - 1st arg
          thrust::make_permutation_iterator( // input - 2nd arg
            tmp_rhod.begin(), 
            tmp_ijk.begin()
          ),
          tmp_real.begin(),                       // output
          arg::_1 * arg::_2 / real_t(rho_stp<real_t>() / si::kilograms * si::cubic_metres)
        ); 
      printf("n (real) po poprawce rhod w init_n\n");
      debug::print(tmp_real);

	// host -> device (includes casting from real_t to uint! and rounding)
	thrust::copy(
          thrust::make_transform_iterator(tmp_real.begin(), arg::_1 + real_t(0.5)),
          thrust::make_transform_iterator(tmp_real.end(), arg::_1 + real_t(0.5)),
          n.begin() + n_part_old
        ); 
      }
      printf("n po skpiowaniu na device i zamianie w int\n");
      debug::print(n);
        
      // detecting possible overflows of n type
      thrust_size_t ix = thrust::max_element(n.begin() + n_part_old, n.end()) - (n.begin() + n_part_old);
      assert(n[ix] < (typename impl::n_t)(-1) / 10000);

    }
  };
};
