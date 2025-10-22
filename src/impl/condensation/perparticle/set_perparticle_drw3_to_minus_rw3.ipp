namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::set_perparticle_drw3_to_minus_rw3(const bool use_stored_rw3)
    {
      thrust_device::vector<real_t> &drw3 = drw3_gp->get();
      if(!use_stored_rw3)
      {
        thrust::transform(
          thrust::make_transform_iterator(rw2.begin(), detail::rw2torw3<real_t>()),
          thrust::make_transform_iterator(rw2.end(), detail::rw2torw3<real_t>()),
          drw3.begin(),
          thrust::negate<real_t>()
        );
      }
      else
      {
        thrust_device::vector<real_t> &rw3 = rw3_gp->get();
        thrust::transform(
          rw3.begin(), rw3.end(),
          drw3.begin(),
          thrust::negate<real_t>()
        );
      }
    }
  };
};