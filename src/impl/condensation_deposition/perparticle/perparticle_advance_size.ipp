namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template<class real_t, bool ice, bool sgs>
      struct RH_hlpr
      {
        std::conditional_t<ice, RH_i<real_t>, RH<real_t> > resolved_RH;

        BOOST_GPU_ENABLED 
        RH_hlpr(RH_formula_t RH_formula):
          resolved_RH(RH_formula)
        {}

        BOOST_GPU_ENABLED 
        real_t operator()(const thrust::tuple<real_t, real_t, real_t, real_t> &tpl) noexcept
        {
          assert(sgs);
          return resolved_RH(thrust::make_tuple(thrust::get<0>(tpl), thrust::get<1>(tpl), thrust::get<2>(tpl))) + thrust::get<3>(tpl);
        }

        BOOST_GPU_ENABLED 
        real_t operator()(const thrust::tuple<real_t, real_t, real_t> &tpl) noexcept
        {
          assert(!sgs);
          return resolved_RH(thrust::make_tuple(thrust::get<0>(tpl), thrust::get<1>(tpl), thrust::get<2>(tpl)));
        }
      };
    };


    template <typename real_t, backend_t device>
    template<bool ice>
    void particles_t<real_t, device>::impl::perparticle_advance_size(
      const real_t &RH_max,
      const bool turb_cond
    ) { 

      namespace arg = thrust::placeholders;

      thrust_device::vector<real_t> &Tp = Tp_gp->get();

      // std::conditional_t<ice, decltype(perparticle_advance_ice_a_c), decltype(perparticle_advance_hlpr<ice>)> advance_func = ice ? advance_ice_a_c<real_t, device, detail::RH_sgs<real_t, true>> : advance_hlpr<ice><real_t, device, detail::RH_sgs<real_t, false>>;

      // advance rw2
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
          perparticle_advance_hlpr<ice>(RH_max, Tp,
            pressure_iter,
            thrust::make_transform_iterator(
              thrust::make_zip_iterator(thrust::make_tuple(
                pressure_iter,
                sstp_tmp_rv.begin(),
                Tp.begin(),
                ssp.begin()
              )),
              detail::RH_hlpr<real_t, ice, true>(opts_init.RH_formula)
            )
          ); 
        else
          perparticle_advance_hlpr<ice>
          (RH_max, Tp,
            pressure_iter,
            thrust::make_transform_iterator(
              thrust::make_zip_iterator(thrust::make_tuple(
                pressure_iter,
                sstp_tmp_rv.begin(),
                Tp.begin()
              )),
              detail::RH_hlpr<real_t, ice, false>(opts_init.RH_formula)
              // detail::RH<real_t>(opts_init.RH_formula)
            )
          ); 
      }
      else
      {
        if(turb_cond)
          perparticle_advance_hlpr<ice>(RH_max, Tp,
            sstp_tmp_p.begin(),
            thrust::make_transform_iterator(
              thrust::make_zip_iterator(thrust::make_tuple(
                sstp_tmp_p.begin(),
                sstp_tmp_rv.begin(),
                Tp.begin(),
                ssp.begin()
              )),  
              detail::RH_hlpr<real_t, ice, true>(opts_init.RH_formula)
            )
          ); 
        else
          perparticle_advance_hlpr<ice>(RH_max, Tp,
            sstp_tmp_p.begin(),
            thrust::make_transform_iterator(
              thrust::make_zip_iterator(thrust::make_tuple(
                sstp_tmp_p.begin(),
                sstp_tmp_rv.begin(),
                Tp.begin()
              )),
              detail::RH_hlpr<real_t, ice, false>(opts_init.RH_formula)
            )
          ); 
      }
    }


    // template <typename real_t, backend_t device>
    // void particles_t<real_t, device>::impl::cond_perparticle_advance_hlpr<ice>(
    //   const real_t &RH_max,
    //   const bool turb_cond
    // ) { 

    //   namespace arg = thrust::placeholders;

    //   thrust_device::vector<real_t> &Tp = Tp_gp->get();

    //   // advance rw2
    //   if(!opts_init.const_p)
    //   {
    //     auto pressure_iter = thrust::make_transform_iterator(
    //       thrust::make_zip_iterator(thrust::make_tuple(
    //         sstp_tmp_rh.begin(),
    //         sstp_tmp_rv.begin(),
    //         Tp.begin()
    //       )),
    //       detail::common__theta_dry__p<real_t>()
    //     );

    //     if(turb_cond)
    //       perparticle_advance_hlpr<ice>(RH_max, Tp,
    //         pressure_iter,
    //         thrust::make_transform_iterator(
    //           thrust::make_zip_iterator(thrust::make_tuple(
    //             pressure_iter,
    //             sstp_tmp_rv.begin(),
    //             Tp.begin(),
    //             ssp.begin()
    //           )),
    //           detail::RH_sgs<real_t>(opts_init.RH_formula)
    //         )
    //       ); 
    //     else
    //       perparticle_advance_hlpr<ice>
    //       (RH_max, Tp,
    //         pressure_iter,
    //         thrust::make_transform_iterator(
    //           thrust::make_zip_iterator(thrust::make_tuple(
    //             pressure_iter,
    //             sstp_tmp_rv.begin(),
    //             Tp.begin()
    //           )),
    //           detail::RH<real_t>(opts_init.RH_formula)
    //         )
    //       ); 
    //   }
    //   else
    //   {
    //     if(turb_cond)
    //       perparticle_advance_hlpr<ice>(RH_max, Tp,
    //         sstp_tmp_p.begin(),
    //         thrust::make_transform_iterator(
    //           thrust::make_zip_iterator(thrust::make_tuple(
    //             sstp_tmp_p.begin(),
    //             sstp_tmp_rv.begin(),
    //             Tp.begin(),
    //             ssp.begin()
    //           )),
    //           detail::RH_sgs<real_t>(opts_init.RH_formula)              
    //         )
    //       ); 
    //     else
    //       perparticle_advance_hlpr<ice>(RH_max, Tp,
    //         sstp_tmp_p.begin(),
    //         thrust::make_transform_iterator(
    //           thrust::make_zip_iterator(thrust::make_tuple(
    //             sstp_tmp_p.begin(),
    //             sstp_tmp_rv.begin(),
    //             Tp.begin()
    //           )),
    //           detail::RH<real_t>(opts_init.RH_formula)              
    //         )
    //       ); 
    //   }
    // }
  };
};
