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

    // init
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_dry(
      const real_t kappa,
      const common::unary_function<real_t> *n_of_lnrd_stp // TODO: kappa-spectrum map
    )
    {
      // memory allocation
      rd3.resize(n_part);
      n.resize(n_part);
      kpa.resize(n_part); 

      // filling kappas
      thrust::fill(kpa.begin(), kpa.end(), kappa);

      // tossing random numbers [0,1] for dry radii
      rand_u01(n_part);

      // values to start the search 
      const real_t rd_min_init = 1e-11, rd_max_init = 1e-3;
      real_t rd_min = rd_min_init, rd_max = rd_max_init;

      // temporary space on the host 
      thrust::host_vector<real_t> tmp(n_part);
      thrust::host_vector<thrust_size_t> tmp_ijk(n_part);
      thrust::host_vector<real_t> &tmp_rhod(tmp_host_real_cell);

      thrust::copy(
	rhod.begin(), rhod.end(), // from
	tmp_rhod.begin()          // to
      );
      thrust::copy(
	ijk.begin(), ijk.end(), // from
	tmp_ijk.begin()         // to
      );

      bool found_optimal_range = 0;

      while (!found_optimal_range)
      {
        namespace arg = thrust::placeholders;

	// rd3 temporarily means logarithm of radius!
	thrust_device::vector<real_t> &lnrd(rd3);

	// shifting from [0,1] to [log(rd_min),log(rd_max)] and storing into rd3
        // each radius randomized only on a small subrange to make the distributions more uniform
        // lnrd is sorted
	thrust::transform(
	  u01.begin(), u01.end(), 
          thrust::make_counting_iterator(real_t(0)),
	  lnrd.begin(), 
	  log(rd_min) + (arg::_2 + arg::_1) * (log(rd_max) - log(rd_min))  / real_t(n_part)
	);

	// filling n with multiplicities
	// (performing it on a local copy as n_of_lnrd_stp may lack __device__ qualifier)
	real_t multiplier = log(rd_max / rd_min) 
          / opts_init.sd_conc_mean
          * opts_init.dx 
          * opts_init.dy 
          * opts_init.dz;

	// device -> host (not needed for omp or cpp ... but happens just once)
	thrust::copy(lnrd.begin(), lnrd.end(), tmp.begin()); 

	// evaluating n_of_lnrd_stp
	thrust::transform(
	  tmp.begin(), tmp.end(), // input 
	  tmp.begin(),            // output
	  detail::eval_and_multiply<real_t>(*n_of_lnrd_stp, multiplier)
	);

        // correcting STP -> actual ambient conditions
        {
          namespace arg = thrust::placeholders;
          using common::earth::rho_stp;

	  thrust::transform(
            tmp.begin(), tmp.end(),            // input - 1st arg
            thrust::make_permutation_iterator( // input - 2nd arg
              tmp_rhod.begin(), 
              tmp_ijk.begin()
            ),
            tmp.begin(),                       // output
            arg::_1 * arg::_2 / real_t(rho_stp<real_t>() / si::kilograms * si::cubic_metres)
          ); 

  	  // host -> device (includes casting from real_t to uint! and rounding)
  	  thrust::copy(
            thrust::make_transform_iterator(tmp.begin(), arg::_1 + real_t(0.5)),
            thrust::make_transform_iterator(tmp.end(), arg::_1 + real_t(0.5)),
            n.begin()); 
        }
        found_optimal_range = 1;

	// chosing an optimal rd_min/rd_max range for a given pdf and grid
	thrust_size_t ix;
	ix = thrust::find_if(n.begin(), n.end(), arg::_1 != 0) - n.begin();
	if (rd_min == rd_min_init) 
	{
          if(ix == n_part)
            std::runtime_error("Initial dry radii distribution outside of the range [1e-11, 1e-3] meters\n");
          if(ix == 0)
            std::runtime_error("Initial dry radii distribution is non-zero for r=1e-11 meters\n");
	  rd_min = exp(lnrd[ix-1]); // adjusting the range
	}
        else if (ix>0)
        {
	  rd_min = exp(lnrd[ix]); // adjusting the range
          found_optimal_range = 0;
        }

	ix = n.rend() - thrust::find_if(n.rbegin(), n.rend(), arg::_1 != 0);
	if (rd_max == rd_max_init) 
	{
          if(ix == n_part)
            std::runtime_error("Initial dry radii distribution is non-zero for r=1e-3 meters\n");
	  rd_max = exp(lnrd[ix+1]); // adjusting the range
          found_optimal_range = 0;
	}
        else if (ix < n_part)
        {
	  rd_max = exp(lnrd[ix]); // adjusting the range
          found_optimal_range = 0;
        }

        if(found_optimal_range)
        {
          // detecting possible overflows of n type
          ix = thrust::max_element(n.begin(), n.end()) - n.begin();
          assert(n[ix] < (typename impl::n_t)(-1) / 10000);

          // converting rd back from logarithms to rd3
          thrust::transform(
            lnrd.begin(),
            lnrd.end(),
            rd3.begin(),
            detail::exp3x<real_t>()
          );
        }
      }
    }
  };
};
