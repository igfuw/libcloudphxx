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
      if (opts_init.sstp_cond == 1) return;
std::cerr << "SAVE" << std::endl;

      const int n = 3;
      thrust_device::vector<real_t>
        *fr[n] = { &rv,          &th,          &rhod        },
        *to[n] = { &sstp_tmp_rv, &sstp_tmp_th, &sstp_tmp_rh };
      for (int ix = 0; ix < n; ++ix)
      {
        thrust::copy(
          fr[ix]->begin(),
          fr[ix]->end(),
          to[ix]->begin()
        ); 
      }
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_sstp()
    {   
std::cerr << "INIT" << std::endl;
      if (opts_init.sstp_cond == 1) return;

      // memory allocation 
      const int n = 3;
      thrust_device::vector<real_t>
        *v[n] = { &sstp_tmp_rv, &sstp_tmp_th, &sstp_tmp_rh };

      for (int ix = 0; ix < n; ++ix)
	v[ix]->resize(n_cell);

      // initialise _old values
      sstp_save();
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::sstp_step(
      const int &step
    )
    {   
std::cerr << "STEP" << std::endl;
      namespace arg = thrust::placeholders;

      const int n = 3;
      thrust_device::vector<real_t>
        *scl[n] = { &rv,          &th,          &rhod        },
        *tmp[n] = { &sstp_tmp_rv, &sstp_tmp_th, &sstp_tmp_rh };

      for (int ix = 0; ix < n; ++ix)
      {
        const int sstp = opts_init.sstp_cond;
	if (step == 0)
	{
	  // sstp_tmp_scl = dscl_adv
	  thrust::transform(
	    scl[ix]->begin(), scl[ix]->end(),     // 1st arg: rv_new
	    tmp[ix]->begin(),                     // 2nd arg: rv_old
	    tmp[ix]->begin(),                     // output (in-place)
	    arg::_1 - arg::_2                     // op: dscl_adv = (scl_new - scl_old)
	  );
	  // scl -= (sstp - 1) * dscl_adv / sstp
	  thrust::transform(
	    scl[ix]->begin(), scl[ix]->end(),     // 1st arg
	    tmp[ix]->begin(),                     // 2nd arg
	    scl[ix]->begin(),                     // output (in-place)
	    arg::_1 - (sstp - 1) * arg::_2 / sstp // op: rv = rv - (sstp - 1) * dscl_adv / sstp
	  );
	}
	else
	{
	  // scl += dscl_adv / sstp
	  thrust::transform(
	    scl[ix]->begin(), scl[ix]->end(),     // 1st arg
	    tmp[ix]->begin(),                     // 2nd arg
	    scl[ix]->begin(),                     // output (in-place)
	    arg::_1 + arg::_2 / sstp              // op: rv = rv + drv_adv / sstp
	  );
	}
      }
    }
  };  
};
