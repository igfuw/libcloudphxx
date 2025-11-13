
namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::apply_noncond_perparticle_sstp_delta()
    {
      namespace arg = thrust::placeholders;

      const int n = 4;

      thrust_device::vector<real_t>
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
          tmp[ix]->begin(), tmp[ix]->end(),
          dlt[ix]->begin(),
          tmp[ix]->begin(),
          arg::_1 + arg::_2 / sstp_cond
        );
      }
    }
  };
};