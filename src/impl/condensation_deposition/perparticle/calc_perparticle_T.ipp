namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::calc_perparticle_T(
    ) { 

      namespace arg = thrust::placeholders;

      thrust_device::vector<real_t> &Tp = Tp_gp->get();

      // calculate perparticle temperature
      if(opts_init.th_dry)
      {
        thrust::transform(
          sstp_tmp_th.begin(), sstp_tmp_th.end(),
          sstp_tmp_rh.begin(),
          Tp.begin(),
          detail::common__theta_dry__T_rhod<real_t>() 
        );  
      }
      else
      {
        thrust::transform(
          sstp_tmp_th.begin(), sstp_tmp_th.end(),
          thrust::make_zip_iterator(thrust::make_tuple(
            sstp_tmp_rv.begin(),
            sstp_tmp_p.begin()
          )),
          Tp.begin(),
          detail::common__theta_std__T_p<real_t>() 
        );
      }
    }
  };
};
