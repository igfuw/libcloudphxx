namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::apply_perparticle_drw3_or_d_ice_mass_to_perparticle_rv_and_th(const bool rw3_changed, const bool ice_mass_changed)
    {
      assert(rw3_changed || ice_mass_changed);

      namespace arg = thrust::placeholders;
      thrust_device::vector<real_t> &Tp = Tp_gp->get();

      if(rw3_changed)
      {
        // get change in liquid water mass
        thrust_device::vector<real_t> &drw3 = drwX_gp->get(); 
        // ice crystals should have drw3==0
        thrust::transform(
          drw3.begin(), drw3.end(),
          drw3.begin(), // in-place!
          arg::_1 * (common::moist_air::rho_w<real_t>() / si::kilograms * si::cubic_metres)
          * real_t(4./3) * pi<real_t>()
        );  

        // calculate rv change from the change in liquid water mass
        thrust::transform(
          drw3.begin(), drw3.end(),
          thrust::make_zip_iterator(thrust::make_tuple(
            sstp_tmp_rh.begin(),
            n.begin(),
            thrust::make_permutation_iterator(dv.begin(), ijk.begin())
          )),
          drw3.begin(),
          detail::massdiff2drv<real_t>(real_t(-1), n_dims)
        );  

        // apply rv change to rv (of this particle only or of all particles, depending on sstp_cond_mix)
        if(opts_init.sstp_cond_mix)
          update_pstate(sstp_tmp_rv, drw3);
        else
          thrust::transform(
            drw3.begin(), drw3.end(),
            sstp_tmp_rv.begin(),
            sstp_tmp_rv.begin(),
            thrust::plus<real_t>()
          );

        // calculate dth from drv
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
      }

      // add the change in ice mass // TODO: very similar to the part for liquid water...
      if(ice_mass_changed)
      {
        thrust_device::vector<real_t> &d_ice_mass = d_ice_mass_percell_gp->get();
        // calculate rv change from the change in solid water mass
        thrust::transform(
          d_ice_mass.begin(), d_ice_mass.end(),
          thrust::make_zip_iterator(thrust::make_tuple(
            sstp_tmp_rh.begin(),
            n.begin(),
            thrust::make_permutation_iterator(dv.begin(), ijk.begin())
          )),
          d_ice_mass.begin(),
          detail::massdiff2drv<real_t>(real_t(-1), n_dims)
        );  

        // apply rv change to rv (of this particle only or of all particles, depending on sstp_cond_mix)
        if(opts_init.sstp_cond_mix)
          update_pstate(sstp_tmp_rv, d_ice_mass);
        else
          thrust::transform(
            d_ice_mass.begin(), d_ice_mass.end(),
            sstp_tmp_rv.begin(),
            sstp_tmp_rv.begin(),
            thrust::plus<real_t>()
          );

        // calculate dth from the change in ice mass
        thrust::transform(
          thrust::make_zip_iterator(thrust::make_tuple(  
            d_ice_mass.begin(),
            Tp.begin(),
            sstp_tmp_th.begin()
          )), 
          thrust::make_zip_iterator(thrust::make_tuple(  
            d_ice_mass.end(),
            Tp.end(),
            sstp_tmp_th.end()
          )), 
          d_ice_mass.begin(),
          detail::dth_dep<real_t>()
        );
      }

      thrust_device::vector<real_t> &dth = rw3_changed ? drwX_gp->get() : d_ice_mass_percell_gp->get();
      if(rw3_changed && ice_mass_changed)
        thrust::transform(
          dth.begin(), dth.end(),
          d_ice_mass_percell_gp->get().begin(),
          dth.begin(),
          thrust::plus<real_t>()
        );

        // apply dth from condensation and deposition
        if(opts_init.sstp_cond_mix)
          update_pstate(sstp_tmp_th, dth);
        else
          thrust::transform(
            dth.begin(), dth.end(),
            sstp_tmp_th.begin(),
            sstp_tmp_th.begin(),
            thrust::plus<real_t>()
          );
    }
  };
};
