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
      // values to start the search 
      const real_t rd_min_init = 1e-11, rd_max_init = 1e-3;
      real_t rd_min = rd_min_init, rd_max = rd_max_init;

      // memory allocation
      rd3.resize(n_part);
      n.resize(n_part);
      kpa.resize(n_part); 

      // filling kappas
      thrust::fill(kpa.begin(), kpa.end(), kappa);

      // tossing random numbers [0,1] for dry radii
      rand_u01(n_part);

      // temporary space on the host 
      thrust::host_vector<real_t> tmp_real(n_part);
      thrust::host_vector<n_t> tmp_n(n_part);
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

        thrust_device::vector<thrust_size_t> &ptr(tmp_device_size_cell);
        thrust::exclusive_scan(count_num.begin(), count_num.end(), ptr.begin()); // number of SDs in cells up to (i-1)

	// shifting from [0,1] to [log(rd_min),log(rd_max)] and storing into rd3
        // each radius randomized only on a small subrange to make the distributions more uniform
        // subranges are specific for each cell
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
/*          thrust::make_zip_iterator(thrust::make_tuple(
            u01.end(), 
            thrust::make_permutation_iterator(count_num.end(), sorted_ijk.end()),     
            thrust::make_counting_iterator(n_part),                                   
            thrust::make_permutation_iterator(ptr.end(), sorted_ijk.end())            
          )),*/
	  lnrd.begin(), 
          calc_lnrd<real_t>(log(rd_min), log(rd_max))
	);

	// filling n with multiplicities
	// (performing it on a local copy as n_of_lnrd_stp may lack __device__ qualifier)
	real_t multiplier = log(rd_max / rd_min) 
          / real_t(opts_init.sd_conc)
          * opts_init.dx 
          * opts_init.dy 
          * opts_init.dz;

	// device -> host (not needed for omp or cpp ... but happens just once)
	thrust::copy(lnrd.begin(), lnrd.end(), tmp_real.begin()); 

	// evaluating n_of_lnrd_stp
	thrust::transform(
	  tmp_real.begin(), tmp_real.end(), // input 
	  tmp_real.begin(),            // output
	  detail::eval_and_oper<real_t>(*n_of_lnrd_stp, multiplier)
	);

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

  	  // host -> host (includes casting from real_t to uint! and rounding)
  	  thrust::copy(
            thrust::make_transform_iterator(tmp_real.begin(), arg::_1 + real_t(0.5)),
            thrust::make_transform_iterator(tmp_real.end(), arg::_1 + real_t(0.5)),
            tmp_n.begin()); 

          // host->device
          thrust::copy(tmp_n.begin(), tmp_n.end(), n.begin());
        }
        found_optimal_range = 1;

	// chosing an optimal rd_min/rd_max range for a given pdf and grid
        // doing it on temp copies of n and lnrd sorted by lnrd
	thrust::copy(lnrd.begin(), lnrd.end(), tmp_real.begin()); 
        thrust::sort_by_key(tmp_real.begin(), tmp_real.end(), tmp_n.begin());

	thrust_size_t ix;
	ix = thrust::find_if(tmp_n.begin(), tmp_n.end(), arg::_1 != 0) - tmp_n.begin();
	if (rd_min == rd_min_init) 
	{
          if(ix == n_part)
            std::runtime_error("Initial dry radii distribution outside of the range [1e-11, 1e-3] meters\n");
          if(ix == 0)
            std::runtime_error("Initial dry radii distribution is non-zero for r=1e-11 meters\n");
	  rd_min = exp(tmp_real[ix-1]); // adjusting the range
	}
        else if (ix>0)
        {
	  rd_min = exp(tmp_real[ix]); // adjusting the range
          found_optimal_range = 0;
        }

	ix = tmp_n.rend() - thrust::find_if(tmp_n.rbegin(), tmp_n.rend(), arg::_1 != 0);
	if (rd_max == rd_max_init) 
	{
          if(ix == n_part)
            std::runtime_error("Initial dry radii distribution is non-zero for r=1e-3 meters\n");
	  rd_max = exp(tmp_real[ix+1]); // adjusting the range
          found_optimal_range = 0;
	}
        else if (ix < n_part)
        {
	  rd_max = exp(tmp_real[ix]); // adjusting the range
          found_optimal_range = 0;
        }

        if(found_optimal_range)
        {
          // detecting possible overflows of n type
          ix = thrust::max_element(tmp_n.begin(), tmp_n.end()) - tmp_n.begin();
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
