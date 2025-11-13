namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    template <class pres_iter, class RH_iter>
    void particles_t<real_t, device>::impl::perparticle_drw2(
      const real_t &RH_max,
      const thrust_device::vector<real_t> &Tp,
      const pres_iter &pi,
      const RH_iter &rhi
    ) { 
      thrust_device::vector<real_t> &lambda_D(lambda_D_gp->get()); 
      thrust_device::vector<real_t> &lambda_K(lambda_K_gp->get()); 
      thrust_device::vector<real_t> &drw2 = drw2_gp->get();

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

      thrust::transform(
        rw2.begin(), rw2.end(),
        thrust::make_zip_iterator(thrust::make_tuple(
          hlpr_zip_iter,
          pi,
          rhi
        )), 
        drw2.begin(),
        detail::advance_rw2<real_t, false>(dt / sstp_cond, RH_max, config.eps_tolerance, config.cond_mlt, config.n_iter)
        );
    }
  };
};