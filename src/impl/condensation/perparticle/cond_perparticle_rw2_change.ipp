namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template<class real_t>
      struct RH_sgs
      {
        RH<real_t> resolved_RH;

        RH_sgs(RH_formula_t RH_formula):
          resolved_RH(RH_formula)
        {}

        BOOST_GPU_ENABLED 
        real_t operator()(const thrust::tuple<real_t, real_t, real_t, real_t> &tpl) 
        {
          return resolved_RH(thrust::make_tuple(thrust::get<0>(tpl), thrust::get<1>(tpl), thrust::get<2>(tpl))) + thrust::get<3>(tpl);
        }
      };
    };


    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::cond_perparticle_rw2_change(
      const real_t &dt,
      const real_t &RH_max,
      const bool turb_cond
    ) { 

      namespace arg = thrust::placeholders;

      thrust_device::vector<real_t> &Tp = Tp_gp->get();

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
          perparticle_advance_rw2(dt, RH_max, Tp,
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
          perparticle_advance_rw2(dt, RH_max, Tp,
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
          perparticle_advance_rw2(dt, RH_max, Tp,
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
          perparticle_advance_rw2(dt, RH_max, Tp,
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