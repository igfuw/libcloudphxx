namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t>
      struct drw2_to_drw3
      {
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &drw2, const real_t &rw2) const
        {
          return real_t(3)/real_t(2) * drw2 * sqrt(rw2);
        }
      };
    }; // namespace detail

    template <typename real_t, backend_t device>
    template <bool use_unconverged_mask>
    void particles_t<real_t, device>::impl::cond_perparticle_drw3_from_drw2()
    {
      thrust_device::vector<real_t> &drw2 = drw2_gp->get(); 
      thrust_device::vector<real_t> &drw3 = drw3_gp->get(); 

      if constexpr (!use_unconverged_mask)
      {
        thrust::transform(
          drw2.begin(), drw2.end(),
          rw2.begin(),
          drw3.begin(),
          detail::drw2_to_drw3<real_t>()
        );
      }
      else
      {
        const auto &unconverged_mask = sstp_cond_unconverged_mask_gp->get();

        thrust::transform_if(
          drw2.begin(), drw2.end(),
          rw2.begin(),
          unconverged_mask.begin(),
          drw3.begin(),
          detail::drw2_to_drw3<real_t>(),
          cuda::std::identity()
        );
      }
    }
  };
};