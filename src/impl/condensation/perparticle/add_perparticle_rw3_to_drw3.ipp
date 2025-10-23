namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    template <int power, bool use_unconverged_mask>
    void particles_t<real_t, device>::impl::add_perparticle_rwX_to_drwX(const bool store_rw3)
    {
      thrust_device::vector<real_t> &drwX = drwX_gp->get();
      if(store_rw3)
      {
        assert(power == 3);
        assert(use_unconverged_mask == false);
        thrust_device::vector<real_t> &rwX = rwX_gp->get();
        thrust::transform(
          rw2.begin(), rw2.end(),
          rwX.begin(),
          detail::rw2torwX<real_t, power>()
        );

        thrust::transform(
          rwX.begin(), rwX.end(),
          drwX.begin(),
          drwX.begin(),
          thrust::plus<real_t>()
        );
      }
      else
      {
        if constexpr (!use_unconverged_mask)
          thrust::transform(
            thrust::make_transform_iterator(rw2.begin(), detail::rw2torwX<real_t, power>()),
            thrust::make_transform_iterator(rw2.end(), detail::rw2torwX<real_t, power>()),
            drwX.begin(),
            drwX.begin(),
            thrust::plus<real_t>()
          );
        else
        {
          const auto &unconverged_mask = sstp_cond_unconverged_mask_gp->get();
          thrust::transform_if(
            thrust::make_transform_iterator(rw2.begin(), detail::rw2torwX<real_t, power>()),
            thrust::make_transform_iterator(rw2.end(), detail::rw2torwX<real_t, power>()),
            drwX.begin(),
            unconverged_mask.begin(),
            drwX.begin(),
            thrust::plus<real_t>(),
            cuda::std::identity()
          );
        }
      }
    }
  };
};