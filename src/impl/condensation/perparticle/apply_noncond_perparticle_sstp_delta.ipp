
namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t>
      
      struct add_dlt_to_tmp
      { 
        BOOST_GPU_ENABLED
        real_t operator() (real_t tmp, thrust::tuple<real_t, int> dlt_sstp_tpl) const noexcept
        {
          return tmp + thrust::get<0>(dlt_sstp_tpl) / thrust::get<1>(dlt_sstp_tpl);
        }
      };
    }

    template <typename real_t, backend_t device>
    template <bool use_unconverged_mask, class it_t>
    void particles_t<real_t, device>::impl::apply_noncond_perparticle_sstp_delta(const it_t sstp_cond_it)
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
        if constexpr (!use_unconverged_mask)
          thrust::transform(
            tmp[ix]->begin(), tmp[ix]->end(),
            thrust::make_zip_iterator(thrust::make_tuple(
              dlt[ix]->begin(),
              sstp_cond_it
            )),
            tmp[ix]->begin(),
            detail::add_dlt_to_tmp<real_t>{}
          );
        else
        {
          const auto &unconverged_mask = sstp_cond_unconverged_mask_gp->get();
          thrust::transform_if(
            tmp[ix]->begin(), tmp[ix]->end(),
            thrust::make_zip_iterator(thrust::make_tuple(
              dlt[ix]->begin(),
              sstp_cond_it
            )),
            // dlt[ix]->begin(),
            unconverged_mask.begin(),
            tmp[ix]->begin(),
            detail::add_dlt_to_tmp<real_t>{},
            cuda::std::identity()
          );
        }
      }
    }
  };
};