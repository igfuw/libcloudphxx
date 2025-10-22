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
      if (!allow_sstp_cond) return;

      const int n = 4;
      thrust_device::vector<real_t>
        *fr[n] = { &rv,          &th,          &rhod,        &p          },
        *to[n] = { &sstp_tmp_rv, &sstp_tmp_th, &sstp_tmp_rh, &sstp_tmp_p };
      for (int ix = 0; ix < ( (opts_init.const_p && opts_init.exact_sstp_cond) ? n : n-1); ++ix) // TODO: var_rho
      {
        if(opts_init.exact_sstp_cond) // per-particle version
          thrust::copy(
            thrust::make_permutation_iterator(fr[ix]->begin(), ijk.begin()),
            thrust::make_permutation_iterator(fr[ix]->begin(), ijk.end()),
            to[ix]->begin()
          ); 
        else // per-cell version
          thrust::copy(
            fr[ix]->begin(),
            fr[ix]->end(),
            to[ix]->begin()
          ); 
      }
    }

    // init _old of new particles, only needed in per-particle sstp cond
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_perparticle_sstp()
    {   
      if (!allow_sstp_cond || !opts_init.exact_sstp_cond) return;

      const int n = 4;
      thrust_device::vector<real_t>
        *fr[n] = { &rv,          &th,          &rhod,        &p          },
        *to[n] = { &sstp_tmp_rv, &sstp_tmp_th, &sstp_tmp_rh, &sstp_tmp_p };
      for (int ix = 0; ix < ( (opts_init.const_p) ? n : n-1); ++ix) // TODO: var_rho
      {
        thrust::copy(
          thrust::make_permutation_iterator(fr[ix]->begin(), ijk.begin()+n_part_old),
          thrust::make_permutation_iterator(fr[ix]->begin(), ijk.end()),
          to[ix]->begin()+n_part_old
        ); 
      }
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::sstp_percell_step(
      const int &step
    )
    {   
      if (sstp_cond == 1) return;

      namespace arg = thrust::placeholders;

      const int n = 3;
      thrust_device::vector<real_t>
        *scl[n] = { &rv,          &th,          &rhod        },
        *tmp[n] = { &sstp_tmp_rv, &sstp_tmp_th, &sstp_tmp_rh };

      for (int ix = 0; ix < (var_rho ? n : n-1); ++ix)
      {
        const real_t sstp = sstp_cond;
        if (step == 0)
        {
          // sstp_tmp_scl = dscl_adv (i.e. delta, i.e. new - old)
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



    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::acquire_arrays_for_perparticle_sstp()
    {        
      reset_guardp(sstp_dlt_rv_gp, tmp_device_real_part); 
      reset_guardp(sstp_dlt_th_gp, tmp_device_real_part); 
      reset_guardp(sstp_dlt_rhod_gp, tmp_device_real_part); 
      if(opts_init.const_p)
        reset_guardp(sstp_dlt_p_gp, tmp_device_real_part); 

      reset_guardp(rw3_gp, tmp_device_real_part);
      reset_guardp(drw3_gp, tmp_device_real_part);
      reset_guardp(Tp_gp, tmp_device_real_part);
    }


    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::release_arrays_for_perparticle_sstp()
    {        
      sstp_dlt_rv_gp.reset();
      sstp_dlt_th_gp.reset();
      sstp_dlt_rhod_gp.reset();
      if(opts_init.const_p)
        sstp_dlt_p_gp.reset();

      rw3_gp.reset();
      drw3_gp.reset();
      Tp_gp.reset();
    }


    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::calculate_noncond_perparticle_sstp_delta()
    { 
      namespace arg = thrust::placeholders;

      const int n = 4;
      thrust_device::vector<real_t>
        *scl[n] = { &rv,          &th,          &rhod,        &p          },
        *tmp[n] = { &sstp_tmp_rv, &sstp_tmp_th, &sstp_tmp_rh, &sstp_tmp_p },
        *dlt[n];

      dlt[0] = &(sstp_dlt_rv_gp->get());
      dlt[1] = &(sstp_dlt_th_gp->get());
      dlt[2] = &(sstp_dlt_rhod_gp->get());
      if(opts_init.const_p)
        dlt[3] = &(sstp_dlt_p_gp->get());

      for (int ix = 0; ix < (opts_init.const_p ? n : n-1); ++ix)
      {
        const real_t sstp = sstp_cond;
        // sstp_tmp_scl = dscl_adv (i.e. delta, i.e. new - old)
        thrust::transform(
          thrust::make_permutation_iterator(scl[ix]->begin(), ijk.begin()), // 1st arg: rv_new
          thrust::make_permutation_iterator(scl[ix]->begin(), ijk.end()), // 1st arg: rv_new
          tmp[ix]->begin(),                     // 2nd arg: rv_old
          dlt[ix]->begin(),                     // output 
          (arg::_1 - arg::_2) / sstp            // op: dscl_adv = (scl_new - scl_old) / sstp
        );
      }
    }


    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::apply_noncond_perparticle_sstp_delta()
    { 
      namespace arg = thrust::placeholders;

      const int n = 4;
      thrust_device::vector<real_t>
        *scl[n] = { &rv,          &th,          &rhod,        &p          },
        *tmp[n] = { &sstp_tmp_rv, &sstp_tmp_th, &sstp_tmp_rh, &sstp_tmp_p },
        *dlt[n];

      dlt[0] = &(sstp_dlt_rv_gp->get());
      dlt[1] = &(sstp_dlt_th_gp->get());
      dlt[2] = &(sstp_dlt_rhod_gp->get());
      if(opts_init.const_p)
        dlt[3] = &(sstp_dlt_p_gp->get());

      for (int ix = 0; ix < (opts_init.const_p ? n : n-1); ++ix)
      {
        thrust::transform(
          tmp[ix]->begin(), tmp[ix]->end(),     // 1st arg
          dlt[ix]->begin(),                     // 2nd arg
          tmp[ix]->begin(),                     // output (in-place)
          arg::_1 + arg::_2                     // op: rv = rv + drv_adv
        );
      }
    }

    // template <typename real_t, backend_t device>
    // void particles_t<real_t, device>::impl::sstp_percell_step_exact(
    //   const int &step
    // )
    // {
    // }


    template <typename real_t, backend_t device>
    // void particles_t<real_t, device>::impl::sstp_percell_step_ssp(
    void particles_t<real_t, device>::impl::apply_perparticle_sgs_supersat(
      const real_t &dt
    )
    {   
      namespace arg = thrust::placeholders;
 
      thrust::transform(
        ssp.begin(), ssp.end(),  
        dot_ssp.begin(),
        ssp.begin(),
        arg::_1 + arg::_2 * dt
      );
    }
  };  
};
