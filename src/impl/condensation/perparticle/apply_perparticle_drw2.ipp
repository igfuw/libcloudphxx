namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    template<bool use_unconverged_mask>
    void particles_t<real_t, device>::impl::apply_perparticle_drw2()
    {
      thrust_device::vector<real_t> &drw2 = drw2_gp->get(); 

      if constexpr (!use_unconverged_mask)
      {
        thrust::transform(
          rw2.begin(), rw2.end(),
          drw2.begin(),
          rw2.begin(),
          thrust::plus<real_t>()
        );
      }
      else
      {
        const auto &unconverged_mask = sstp_cond_unconverged_mask_gp->get();

        thrust::transform_if(
          drw2.begin(), drw2.end(),
          rw2.begin(),
          unconverged_mask.begin(),
          rw2.begin(),
          thrust::plus<real_t>(),
          cuda::std::identity()
        );
      }
    }
  };
};