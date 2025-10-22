namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::add_perparticle_rw3_to_drw3(const bool store_rw3)
    {
      thrust_device::vector<real_t> &drw3 = drw3_gp->get();
      if(store_rw3)
      {
        thrust_device::vector<real_t> &rw3 = rw3_gp->get();
        thrust::transform(
          rw2.begin(), rw2.end(),
          rw3.begin(),
          detail::rw2torw3<real_t>()
        );

        thrust::transform(
          rw3.begin(), rw3.end(),
          drw3.begin(),
          drw3.begin(),
          thrust::plus<real_t>()
        );
      }
      else
      {
        thrust::transform(
          thrust::make_transform_iterator(rw2.begin(), detail::rw2torw3<real_t>()),
          thrust::make_transform_iterator(rw2.end(), detail::rw2torw3<real_t>()),
          drw3.begin(),
          drw3.begin(),
          thrust::plus<real_t>()
        );
      }
    }
  };
};