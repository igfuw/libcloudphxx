namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::set_perparticle_d_ice_mass_to_minus_ice_mass(const bool use_stored_ice_mass)
    {
      thrust_device::vector<real_t> &d_ice_mass = d_ice_mass_gp->get();
      if(!use_stored_ice_mass)
      {
        auto ice_mass_it = 
          thrust::make_transform_iterator(
            thrust::make_zip_iterator(thrust::make_tuple(ice_a.begin(), ice_c.begin(), ice_rho.begin())),
            detail::ice_mass<real_t>()
          );

        thrust::transform(
          ice_mass_it, ice_mass_it + n_part,
          d_ice_mass.begin(),
          thrust::negate<real_t>()
        );
      }
      else
      {
        thrust_device::vector<real_t> &ice_mass = ice_mass_gp->get();
        thrust::transform(
          ice_mass.begin(), ice_mass.end(),
          d_ice_mass.begin(),
          thrust::negate<real_t>()
        );
      }
    }
  };
};

