namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    template <bool ice, class pres_iter, class RH_iter>
    void particles_t<real_t, device>::impl::perparticle_advance_hlpr(
      const real_t &RH_max,
      const thrust_device::vector<real_t> &Tp,
      const pres_iter &pi,
      const RH_iter &rhi
    ) { 
      thrust_device::vector<real_t> &lambda_D(lambda_D_gp->get()); 
      thrust_device::vector<real_t> &lambda_K(lambda_K_gp->get()); 

      namespace arg = thrust::placeholders;

      auto hlpr_zip_iter = thrust::make_zip_iterator(thrust::make_tuple(
        sstp_tmp_rh.begin(),
        sstp_tmp_rv.begin(),
        Tp.begin(),
        thrust::make_transform_iterator(
          Tp.begin(),
          detail::common__vterm__visc<real_t>()
        ),
        rd3.begin(),
        kpa.begin(),
        vt.begin(),
        thrust::make_permutation_iterator(lambda_D.begin(), ijk.begin()),
        thrust::make_permutation_iterator(lambda_K.begin(), ijk.begin())
      ));

      if constexpr(!ice)
      {
        thrust::transform(
          rw2.begin(), rw2.end(),
          thrust::make_zip_iterator(thrust::make_tuple(
            hlpr_zip_iter,
            pi,
            rhi
          )), 
          rw2.begin(),
          detail::advance_rw2<real_t>(dt / sstp_cond, RH_max, config.eps_tolerance, config.cond_mlt, config.n_iter)
        );
      }
      else // ice
      {
        auto hlpr_zip_iter_ice = thrust::make_zip_iterator(thrust::make_tuple(
          ice_a.begin(),
          ice_c.begin()
        ));

        thrust::transform(
          hlpr_zip_iter_ice,
          hlpr_zip_iter_ice + n_part,
          thrust::make_zip_iterator(thrust::make_tuple(
            hlpr_zip_iter,
            pi,
            rhi
          )), 
          hlpr_zip_iter_ice,
          detail::advance_ice_ac<real_t>(dt / sstp_cond, RH_max)
        );
      }    
    }
  };
};
