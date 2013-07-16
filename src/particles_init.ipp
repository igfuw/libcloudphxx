// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief initialisation routine for super droplets
  */

#include <iostream>
#include <algorithm>

#include "detail/urand.hpp"
#include "detail/thrust.hpp"
#include "detail/functors_host.hpp"
#include "detail/functors_device.hpp"

#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/extrema.h>

namespace libcloudphxx
{
  namespace lgrngn
  {
    // init
    template <typename real_t, int device>
    void particles<real_t, device>::init(
      unary_function<real_t> *n_of_lnrd
    )
    {
      using namespace thrust::placeholders;

      // tossing random numbers [0,1]
      pimpl->urand(pimpl->n_part);

      // sorting them (does not harm and makes rd_min/rd_max search simpler)
      thrust::sort(pimpl->u01.begin(), pimpl->u01.end());

      // values to start the search 
      const real_t rd_min_init = 1e-10, rd_max_init = 1e-5;
      real_t rd_min = rd_min_init, rd_max = rd_max_init;

      // temporary space on the host 
      thrust::host_vector<real_t> tmp(pimpl->n_part);

      while (true)
      {
	// pimpl->rd3 temporarily means logarithm of radius!
	thrust_device::vector<real_t> &lnrd(pimpl->rd3);

	// shifting from [0,1] to [log(rd_min),log(rd_max)] and storing into rd3
	thrust::transform(
	  pimpl->u01.begin(), 
	  pimpl->u01.end(), 
	  lnrd.begin(), 
	  log(rd_min) + _1 * (log(rd_max) - log(rd_min)) 
	);
 
	// filling pimpl->n with multiplicities
	// (performing it on a local copy as n_of_lnrd may lack __device__ qualifier)
// TODO
real_t dx = 1, dy = 1, dz = 1;
	real_t multiplier = log(rd_max / rd_min) / pimpl->sd_conc_mean * dx * dy * dz;

	// device -> host (not needed for omp or cpp ... but happens just once)
	thrust::copy(lnrd.begin(), lnrd.end(), tmp.begin()); 
	// evaluating n_of_lnrd
	thrust::transform(
	  tmp.begin(),
	  tmp.end(),
	  tmp.begin(),
	  detail::eval_and_multiply<real_t>(*n_of_lnrd, multiplier)
	);
	// host -> device (includes casting from real_t to uint!)
	thrust::copy(tmp.begin(), tmp.end(), pimpl->n.begin()); 

	// chosing an optimal rd_min/rd_max range for a given pdf and grid
	thrust_size_t ix;
	ix = thrust::find_if(pimpl->n.begin(), pimpl->n.end(), _1 != 0) - pimpl->n.begin();
	if (rd_min == rd_min_init) 
	{
	  assert(ix != pimpl->n_part); // zeros everywhere
	  assert(ix != 0); // rd_min was not small enough for pdf(rd_min) to be zero
	  rd_min = exp(lnrd[ix - 1]); // adjusting the range
	}
	ix = pimpl->n.rend() - thrust::find_if(pimpl->n.rbegin(), pimpl->n.rend(), _1 != 0);
	if (rd_max == rd_max_init) 
	{
	  assert(ix != pimpl->n_part); // rd_max was not big enough for pdf(rd_max) to be zero
	  rd_max = exp(lnrd[ix + 1]); // adjusting the range
	  continue;
	}

	// detecting possible overflows of pimpl->n type
	ix = thrust::max_element(pimpl->n.begin(), pimpl->n.end()) - pimpl->n.begin();
	assert(pimpl->n[ix] < (typename impl::n_t)(-1) / 10000);

	// converting rd back from logarithms to rd3
	thrust::transform(
	  lnrd.begin(),
	  lnrd.end(),
	  pimpl->rd3.begin(),
	  detail::exp3x<real_t>()
	);

	break;
      }
    }
  };
};
