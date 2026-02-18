namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::depo_perparticle_advance_ice_a_c(
      const real_t &RH_max,
      const bool turb_cond
    ) { 

      // TODO: it works with turb cond, but sanity checks think it doesnt?
      // NOTE: also the same as advance rw2, unify them!
      namespace arg = thrust::placeholders;

      thrust_device::vector<real_t> &Tp = Tp_gp->get();

      if(!opts_init.const_p)
      {
        auto pressure_iter = thrust::make_transform_iterator(
          thrust::make_zip_iterator(thrust::make_tuple(
            sstp_tmp_rh.begin(),
            sstp_tmp_rv.begin(),
            Tp.begin()
          )),
          detail::common__theta_dry__p<real_t>()
        );

        if(turb_cond)
          perparticle_advance_ice_a_c(RH_max, Tp,
            pressure_iter,
            thrust::make_transform_iterator(
              thrust::make_zip_iterator(thrust::make_tuple(
                pressure_iter,
                sstp_tmp_rv.begin(),
                Tp.begin(),
                ssp.begin()
              )),
              detail::RH_sgs<real_t>(opts_init.RH_formula)
            )
          ); 
        else
          perparticle_advance_ice_a_c
          (RH_max, Tp,
            pressure_iter,
            thrust::make_transform_iterator(
              thrust::make_zip_iterator(thrust::make_tuple(
                pressure_iter,
                sstp_tmp_rv.begin(),
                Tp.begin()
              )),
              detail::RH<real_t>(opts_init.RH_formula)
            )
          ); 
      }
      else
      {
        if(turb_cond)
          perparticle_advance_ice_a_c(RH_max, Tp,
            sstp_tmp_p.begin(),
            thrust::make_transform_iterator(
              thrust::make_zip_iterator(thrust::make_tuple(
                sstp_tmp_p.begin(),
                sstp_tmp_rv.begin(),
                Tp.begin(),
                ssp.begin()
              )),
              detail::RH_sgs<real_t>(opts_init.RH_formula)              
            )
          ); 
        else
          perparticle_advance_ice_a_c(RH_max, Tp,
            sstp_tmp_p.begin(),
            thrust::make_transform_iterator(
              thrust::make_zip_iterator(thrust::make_tuple(
                sstp_tmp_p.begin(),
                sstp_tmp_rv.begin(),
                Tp.begin()
              )),
              detail::RH<real_t>(opts_init.RH_formula)              
            )
          ); 
      }
    }
  };
};
