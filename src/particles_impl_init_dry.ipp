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
      /// @brief returns real_t(exp(3*x))
      template <typename real_t>
      struct exp3x
      { 
	BOOST_GPU_ENABLED 
	real_t operator()(real_t x) 
	{ 
	  return exp(3*x); 
	} 
      };

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

        real_t operator()(real_t x)  
        {
          return mul * fun(x); 
        }
      };  
    };
    
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

    // init
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_dry(
      const real_t kappa,
      const common::unary_function<real_t> *n_of_lnrd_stp // TODO: kappa-spectrum map
    )
    {
      // probing the spectrum to find rd_min-rd_max range
      // values to start the search 
      const real_t rd_min_init = 1e-11, rd_max_init = 1e-3;
      real_t rd_min = rd_min_init, rd_max = rd_max_init;

      bool found_optimal_range = false;
      real_t multiplier;
      while (!found_optimal_range)
      {
	multiplier = log(rd_max / rd_min) 
	  / opts_init.sd_conc
	  * (n_dims == 0
	    ? dv[0]
	    : (opts_init.dx * opts_init.dy * opts_init.dz)
	  );
        impl::n_t 
          n_min = (*n_of_lnrd_stp)(log(rd_min)) * multiplier, 
          n_max = (*n_of_lnrd_stp)(log(rd_max)) * multiplier;

        if (rd_min == rd_min_init && n_min != 0) 
          throw std::runtime_error("Initial dry radii distribution is non-zero for rd_min_init");
        if (rd_max == rd_max_init && n_max != 0) 
          throw std::runtime_error("Initial dry radii distribution is non-zero for rd_max_init");
        
        if      (n_min == 0) rd_min *= 1.1;
        else if (n_max == 0) rd_max /= 1.1;
        else found_optimal_range = true;
      }

      // filling kappas
      thrust::fill(kpa.begin(), kpa.end(), kappa);

      // tossing random numbers [0,1] for dry radii
      rand_u01(n_part);

      // temporary space on the host 
      thrust::host_vector<real_t> tmp_real(n_part);
      thrust::host_vector<thrust_size_t> tmp_ijk(n_part);
      thrust::host_vector<real_t> &tmp_rhod(tmp_host_real_cell);

      thrust::copy(
	rhod.begin(), rhod.end(), // from
	tmp_rhod.begin()          // to
      );
      // first n_part SDs are to init
      thrust::copy(
	ijk.begin(), ijk.begin() + n_part, // from
	tmp_ijk.begin()         // to
      );

      // rd3 temporarily means logarithm of radius!
      thrust_device::vector<real_t> &lnrd(rd3);
      
      // ! count_num and sorted_ijk were filled in init_xyz !
      thrust_device::vector<thrust_size_t> &ptr(tmp_device_size_cell);
      thrust::exclusive_scan(count_num.begin(), count_num.end(), ptr.begin()); // number of SDs in cells up to (i-1)
      
      // shifting from [0,1] to [log(rd_min),log(rd_max)] and storing into rd3
      // each log(radius) randomized only on a small subrange to make the distributions more uniform
      // particles are sorted by cell number (see particles_impl_init_xyz), uniform distribution in each cell
      // lnrd is not sorted
      thrust::transform(
        thrust::make_zip_iterator(thrust::make_tuple(
          u01.begin(),                                                              // random number
          thrust::make_permutation_iterator(count_num.begin(), sorted_ijk.begin()), // number of SDs in the cell
          thrust::make_counting_iterator(0),                                        // sequence to iterate over distribution
          thrust::make_permutation_iterator(ptr.begin(), sorted_ijk.begin())        // number of SDs in cells up to this one
        )),
        thrust::make_zip_iterator(thrust::make_tuple(
          u01.begin(),                                                              // random number
          thrust::make_permutation_iterator(count_num.begin(), sorted_ijk.begin()), // number of SDs in the cell
          thrust::make_counting_iterator(0),                                        // sequence to iterate over distribution
          thrust::make_permutation_iterator(ptr.begin(), sorted_ijk.begin())        // number of SDs in cells up to this one
        )) + n_part,
        lnrd.begin(), 
        calc_lnrd<real_t>(log(rd_min), log(rd_max))
      );
      
      // filling n with multiplicities
      // (performing it on a local copy as n_of_lnrd_stp may lack __device__ qualifier)
      // device -> host (not needed for omp or cpp ... but happens just once)
      thrust::copy(lnrd.begin(), lnrd.begin() + n_part, tmp_real.begin()); 
      
      // evaluating n_of_lnrd_stp
      thrust::transform(
        tmp_real.begin(), tmp_real.begin() + n_part, // input 
        tmp_real.begin(),                 // output
        detail::eval_and_multiply<real_t>(*n_of_lnrd_stp, multiplier)
      );

      // correcting STP -> actual ambient conditions
      {
        namespace arg = thrust::placeholders;
        using common::earth::rho_stp;

        thrust::transform(
          tmp_real.begin(), tmp_real.begin() + n_part,            // input - 1st arg
          thrust::make_permutation_iterator( // input - 2nd arg
            tmp_rhod.begin(), 
            tmp_ijk.begin()
          ),
          tmp_real.begin(),                       // output
          arg::_1 * arg::_2 / real_t(rho_stp<real_t>() / si::kilograms * si::cubic_metres)
        ); 

	// host -> device (includes casting from real_t to uint! and rounding)
	thrust::copy(
          thrust::make_transform_iterator(tmp_real.begin(), arg::_1 + real_t(0.5)),
          thrust::make_transform_iterator(tmp_real.begin() + n_part, arg::_1 + real_t(0.5)),
          n.begin()
        ); 
      }
        
      // detecting possible overflows of n type
      thrust_size_t ix = thrust::max_element(n.begin(), n.begin() + n_part) - n.begin();
      assert(n[ix] < (typename impl::n_t)(-1) / 10000);

      // converting rd back from logarithms to rd3
      thrust::transform(
        lnrd.begin(),
        lnrd.begin() + n_part,
        rd3.begin(),
        detail::exp3x<real_t>()
      );
    }
  };
};
