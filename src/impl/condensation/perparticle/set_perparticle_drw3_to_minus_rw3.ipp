namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    template <int power, bool use_unconverged_mask>
    void particles_t<real_t, device>::impl::set_perparticle_drwX_to_minus_rwX(const bool use_stored_rw3)
    {
      thrust_device::vector<real_t> &drwX = drwX_gp->get();
      if(!use_stored_rw3)
      {
        if(!use_unconverged_mask)
          thrust::transform(
            thrust::make_transform_iterator(rw2.begin(), detail::rw2torwX<real_t, power>()),
            thrust::make_transform_iterator(rw2.end(), detail::rw2torwX<real_t, power>()),
            drwX.begin(),
            thrust::negate<real_t>()
          );
        else
        {
          auto &unconverged_mask = sstp_cond_unconverged_mask_gp->get();
          thrust::transform_if(
            thrust::make_transform_iterator(rw2.begin(), detail::rw2torwX<real_t, power>()),
            thrust::make_transform_iterator(rw2.end(), detail::rw2torwX<real_t, power>()),
            drwX.begin(),
            unconverged_mask.begin(),
            thrust::negate<real_t>(),
            thrust::identity<bool>()
          );
        }
      }
      else
      {
        assert(power == 3);
        assert(use_converged_mask == false);
        thrust_device::vector<real_t> &rw3 = rwX_gp->get();
        thrust::transform(
          rw3.begin(), rw3.end(),
          drwX.begin(),
          thrust::negate<real_t>()
        );
      }
    }
  };
};