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

        // add the change in ice mass
        if(ice_mass_changed)
          thrust::transform(
            drw3.begin(), drw3.end(),
            d_ice_mass_gp->get().begin(),
            drw3.begin(),
            thrust::plus<real_t>()
          );
      }
      // else, theres no liquid water mass change, hence d_ice_mass already has ice mass change

      // calculate rv change from the change in liquid and solid water mass
      thrust_device::vector<real_t> &drv = rw3_changed ? drwX_gp->get() : d_ice_mass_gp->get();  // NOTE: this leads to in-place changes of these vectors.
      thrust::transform(
        drv.begin(), drv.end(),
        thrust::make_zip_iterator(thrust::make_tuple(
          sstp_tmp_rh.begin(),
          n.begin(),
          thrust::make_permutation_iterator(dv.begin(), ijk.begin())
        )),
        drv.begin(),
        detail::rw3diff2drv<real_t>(real_t(-1), n_dims)
      );  

      // apply rv change to rv (of this particle only or of all particles, depending on sstp_cond_mix)
      if(opts_init.sstp_cond_mix)
        update_pstate(sstp_tmp_rv, drv);
      else
        thrust::transform(
          drv.begin(), drv.end(),
          sstp_tmp_rv.begin(),
          sstp_tmp_rv.begin(),
          thrust::plus<real_t>()
        );

      // calculate dth from drv and apply it
      thrust::transform(
        thrust::make_zip_iterator(thrust::make_tuple(  
          drv.begin(),
          Tp.begin(),
          sstp_tmp_th.begin()
        )), 
        thrust::make_zip_iterator(thrust::make_tuple(  
          drv.end(),
          Tp.end(),
          sstp_tmp_th.end()
        )), 
        drv.begin(),
        detail::dth<real_t>()
      );

      if(opts_init.sstp_cond_mix)
        update_pstate(sstp_tmp_th, drv);
      else
        thrust::transform(
          drv.begin(), drv.end(),
          sstp_tmp_th.begin(),
          sstp_tmp_th.begin(),
          thrust::plus<real_t>()
        );
    }
  };
};
