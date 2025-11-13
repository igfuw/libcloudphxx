namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::apply_perparticle_drw3_to_perparticle_rv_and_th()
    {
      namespace arg = thrust::placeholders;
      thrust_device::vector<real_t> &drw3 = drw3_gp->get(); 
      thrust_device::vector<real_t> &Tp = Tp_gp->get();

      thrust::transform(
        drw3.begin(), drw3.end(),
        thrust::make_zip_iterator(thrust::make_tuple(
          sstp_tmp_rh.begin(),
          n.begin(),
          thrust::make_permutation_iterator(dv.begin(), ijk.begin())
        )),
        drw3.begin(),
        detail::rw3diff2drv<real_t>(
          - common::moist_air::rho_w<real_t>() / si::kilograms * si::cubic_metres
          * real_t(4./3) * pi<real_t>(), n_dims
        )
      );  

      if(opts_init.sstp_cond_mix)
        update_pstate(sstp_tmp_rv, drw3);
      else
        thrust::transform(
          drw3.begin(), drw3.end(),
          sstp_tmp_rv.begin(),
          sstp_tmp_rv.begin(),
          thrust::plus<real_t>()
        );

      thrust::transform(
        thrust::make_zip_iterator(thrust::make_tuple(  
          drw3.begin(),
          Tp.begin(),
          sstp_tmp_th.begin()
        )), 
        thrust::make_zip_iterator(thrust::make_tuple(  
          drw3.end(),
          Tp.end(),
          sstp_tmp_th.end()
        )), 
        drw3.begin(),
        detail::dth<real_t>()
      );

      if(opts_init.sstp_cond_mix)
        update_pstate(sstp_tmp_th, drw3);
      else
        thrust::transform(
          drw3.begin(), drw3.end(),
          sstp_tmp_th.begin(),
          sstp_tmp_th.begin(),
          thrust::plus<real_t>()
        );
    }
  };
};