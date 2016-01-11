// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <thrust/iterator/transform_iterator.h>

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::sstp_save()
    {
      const int n = 3;
      thrust_device::vector<real_t>
        *fr[n] = { &rv,          &th,          &rhod        },
        *to[n] = { &sstp_tmp_rv, &sstp_tmp_th, &sstp_tmp_rh };
      for (int ix = 0; ix < n; ++ix) // TODO: var_rho
      {
        thrust::copy(
          thrust::make_permutation_iterator(fr[ix]->begin(), ijk.begin()),
          thrust::make_permutation_iterator(fr[ix]->begin(), ijk.end()),
          to[ix]->begin()
        ); 
      }
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_sstp()
    {   
      // initialise _old values
      sstp_save();
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::sstp_step(
      const int &step,
      const bool &var_rho // if rho varied and need to be updated
    )
    {   
      if (opts_init.sstp_cond == 1) return;

      namespace arg = thrust::placeholders;

      const int n = 3;
      thrust_device::vector<real_t>
        *scl[n] = { &rv,          &th,          &rhod        },
        *tmp[n] = { &sstp_tmp_rv, &sstp_tmp_th, &sstp_tmp_rh },
        *dlt[n] = { &tmp_device_real_part, &tmp_device_real_part1, &tmp_device_real_part2 };

      for (int ix = 0; ix < (var_rho ? n : n-1); ++ix)
      {
        const real_t sstp = opts_init.sstp_cond;
      	if (step == 0)
      	{
      	  // sstp_tmp_scl = dscl_adv (i.e. delta, i.e. new - old)
      	  thrust::transform(
      	    thrust::make_permutation_iterator(scl[ix]->begin(), ijk.begin()), // 1st arg: rv_new
      	    thrust::make_permutation_iterator(scl[ix]->begin(), ijk.end()), // 1st arg: rv_new
      	    tmp[ix]->begin(),                     // 2nd arg: rv_old
      	    dlt[ix]->begin(),                     // output 
      	    (arg::_1 - arg::_2) / sstp            // op: dscl_adv = (scl_new - scl_old) / sstp
      	  );
      	  // scl -= (sstp - 1) * dscl_adv / sstp
      	  thrust::transform(
      	    thrust::make_permutation_iterator(scl[ix]->begin(), ijk.begin()), // 1st arg: rv_new
      	    thrust::make_permutation_iterator(scl[ix]->begin(), ijk.end()), // 1st arg: rv_new
      	    dlt[ix]->begin(),                     // 2nd arg
      	    tmp[ix]->begin(),                     // output 
      	    arg::_1 - (sstp - 1) * arg::_2        // op: rv = rv - (sstp - 1) * dscl_adv
      	  );
      	}
      	else
      	{
      	  // scl += dscl_adv / sstp
      	  thrust::transform(
      	    tmp[ix]->begin(), tmp[ix]->end(),     // 1st arg
      	    dlt[ix]->begin(),                     // 2nd arg
      	    tmp[ix]->begin(),                     // output (in-place)
      	    arg::_1 + arg::_2                     // op: rv = rv + drv_adv
      	  );
      	}
      }
    }
  };  
};
