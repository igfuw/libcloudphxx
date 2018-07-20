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
    void particles_t<real_t, device>::impl::sstp_save_chem()
    {
      if (opts_init.sstp_chem == 1) return;

      const int n = 6;
      thrust_device::vector<real_t>
        *fr[n] = { &ambient_chem[(chem_species_t)0], &ambient_chem[(chem_species_t)1], &ambient_chem[(chem_species_t)2],
                   &ambient_chem[(chem_species_t)3], &ambient_chem[(chem_species_t)4], &ambient_chem[(chem_species_t)5]},

        *to[n] = { &sstp_tmp_chem_0, &sstp_tmp_chem_1, &sstp_tmp_chem_2, &sstp_tmp_chem_3, &sstp_tmp_chem_4, &sstp_tmp_chem_5};
      for (int ix = 0; ix < n; ++ix) // TODO: var_rho
      {
        thrust::copy(
          fr[ix]->begin(),
          fr[ix]->end(),
          to[ix]->begin()
        ); 
      }
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_sstp_chem()
    {   
      if (opts_init.sstp_chem == 1) return;

      // memory allocation 
      sstp_tmp_chem_0.resize(n_cell);
      sstp_tmp_chem_1.resize(n_cell);
      sstp_tmp_chem_2.resize(n_cell);
      sstp_tmp_chem_3.resize(n_cell);
      sstp_tmp_chem_4.resize(n_cell);
      sstp_tmp_chem_5.resize(n_cell);

      // initialise _old values
      sstp_save_chem();
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::sstp_step_chem(
      const int &step
    )
    {   
      if (opts_init.sstp_chem == 1) return;

      namespace arg = thrust::placeholders;

      const int n = 6;

      thrust_device::vector<real_t>
        *scl[n] = { &ambient_chem[(chem_species_t)0], &ambient_chem[(chem_species_t)1], &ambient_chem[(chem_species_t)2],
                    &ambient_chem[(chem_species_t)3], &ambient_chem[(chem_species_t)4], &ambient_chem[(chem_species_t)5]},

        *tmp[n] = { &sstp_tmp_chem_0, &sstp_tmp_chem_1, &sstp_tmp_chem_2, &sstp_tmp_chem_3, &sstp_tmp_chem_4, &sstp_tmp_chem_5};

      for (int ix = 0; ix < (var_rho ? n : n-1); ++ix)
      {
        const real_t sstp = opts_init.sstp_chem;
	if (step == 0)
	{
	  // sstp_tmp_scl = dscl_adv (i.e. delta, i.e. new - old)
	  thrust::transform(
	    scl[ix]->begin(), scl[ix]->end(),     // 1st arg: gas_new
	    tmp[ix]->begin(),                     // 2nd arg: gas_old
	    tmp[ix]->begin(),                     // output (in-place)
	    arg::_1 - arg::_2                     // op: dscl_adv = (scl_new - scl_old)
	  );

	  // scl -= (sstp - 1) * dscl_adv / sstp
	  thrust::transform(
	    scl[ix]->begin(), scl[ix]->end(),     // 1st arg
	    tmp[ix]->begin(),                     // 2nd arg
	    scl[ix]->begin(),                     // output (in-place)
	    arg::_1 - (sstp - 1) * arg::_2 / sstp // op: gas = gas - (sstp - 1) * dscl_adv / sstp
	  );

	}
	else
	{
	  // scl += dscl_adv / sstp
	  thrust::transform(
	    scl[ix]->begin(), scl[ix]->end(),     // 1st arg
	    tmp[ix]->begin(),                     // 2nd arg
	    scl[ix]->begin(),                     // output (in-place)
	    arg::_1 + arg::_2 / sstp              // op: gas = gas + drv_adv / sstp
	  );
	}
      }
    }
  };  
};
