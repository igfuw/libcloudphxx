
namespace libcloudphxx
{
  namespace lgrngn
  {
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
        thrust::transform(
          thrust::make_permutation_iterator(scl[ix]->begin(), ijk.begin()),
          thrust::make_permutation_iterator(scl[ix]->begin(), ijk.end()),
          tmp[ix]->begin(),
          dlt[ix]->begin(),
          (arg::_1 - arg::_2) // / sstp
        );
    }
    }
  };
};