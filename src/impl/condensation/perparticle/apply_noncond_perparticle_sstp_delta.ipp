
namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    template <bool use_unconverged_mask>
    void particles_t<real_t, device>::impl::apply_noncond_perparticle_sstp_delta(const real_t &multiplier)
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
        if constexpr (!use_unconverged_mask)
          thrust::transform(
            tmp[ix]->begin(), tmp[ix]->end(),
            dlt[ix]->begin(),
            tmp[ix]->begin(),
            arg::_1 + multiplier * arg::_2
          );
        else
        {
          const auto &unconverged_mask = sstp_cond_unconverged_mask_gp->get();
          thrust::transform_if(
            tmp[ix]->begin(), tmp[ix]->end(),
            dlt[ix]->begin(),
            unconverged_mask.begin(),
            tmp[ix]->begin(),
            arg::_1 + multiplier * arg::_2,
            cuda::std::identity()
          );
        }
      }
    }
  };
};