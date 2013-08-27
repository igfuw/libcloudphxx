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
	__device__ 
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
    template <typename real_t, int device>
    void particles<real_t, device>::impl::init_dry(
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

      // sorting them (does not harm and makes rd_min/rd_max search simpler)
      thrust::sort(u01.begin(), u01.end());

      // values to start the search 
      const real_t rd_min_init = 1e-10, rd_max_init = 1e-5;
      real_t rd_min = rd_min_init, rd_max = rd_max_init;

      // temporary space on the host 
      thrust::host_vector<real_t> tmp(n_part);

      while (true)
      {
        using namespace thrust::placeholders;

	// rd3 temporarily means logarithm of radius!
	thrust_device::vector<real_t> &lnrd(rd3);

	// shifting from [0,1] to [log(rd_min),log(rd_max)] and storing into rd3
	thrust::transform(
	  u01.begin(), 
	  u01.end(), 
	  lnrd.begin(), 
	  log(rd_min) + _1 * (log(rd_max) - log(rd_min)) 
	);
 
	// filling n with multiplicities
	// (performing it on a local copy as n_of_lnrd_stp may lack __device__ qualifier)
	real_t multiplier = log(rd_max / rd_min) 
          / opts.sd_conc_mean 
          * opts.dx 
          * opts.dy 
          * opts.dz;

	// device -> host (not needed for omp or cpp ... but happens just once)
	thrust::copy(lnrd.begin(), lnrd.end(), tmp.begin()); 
	// evaluating n_of_lnrd_stp
	thrust::transform(
	  tmp.begin(),
	  tmp.end(),
	  tmp.begin(),
	  detail::eval_and_multiply<real_t>(*n_of_lnrd_stp, multiplier)
	);
        // correcting STP -> actual conditions
        {
          thrust::host_vector<real_t> &rhod_host(tmp_host_real_cell);

	  thrust::copy(
	    rhod.begin(), rhod.end(), // from
	    rhod_host.begin()         // to
	  );

          using namespace thrust::placeholders;
          using common::earth::rho_stp;

	  thrust::transform(
            tmp.begin(), tmp.end(),  // input - 1st arg
            rhod_host.begin(),       // input - 2nd arg
            tmp.begin(),             // output
            _1 * _2 / (rho_stp<real_t>() / si::kilograms * si::cubic_metres)
          ); 
        }
	// host -> device (includes casting from real_t to uint!)
	thrust::copy(tmp.begin(), tmp.end(), n.begin()); 

	// chosing an optimal rd_min/rd_max range for a given pdf and grid
	thrust_size_t ix;
	ix = thrust::find_if(n.begin(), n.end(), _1 != 0) - n.begin();
	if (rd_min == rd_min_init) 
	{
	  assert(ix != n_part); // zeros everywhere
	  assert(ix != 0); // rd_min was not small enough for pdf(rd_min) to be zero
	  rd_min = exp(lnrd[ix - 1]); // adjusting the range
	}
	ix = n.rend() - thrust::find_if(n.rbegin(), n.rend(), _1 != 0);
	if (rd_max == rd_max_init) 
	{
	  assert(ix != n_part); // rd_max was not big enough for pdf(rd_max) to be zero
	  rd_max = exp(lnrd[ix + 1]); // adjusting the range
	  continue;
	}

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

	break;
      }
    }
  };
};
