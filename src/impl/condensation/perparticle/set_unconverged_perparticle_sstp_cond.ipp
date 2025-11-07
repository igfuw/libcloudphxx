namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::set_unconverged_perparticle_sstp_cond(const unsigned int &n) noexcept
    {
      auto &perparticle_sstp_cond = perparticle_sstp_cond_gp->get();
      const auto &unconverged_mask = sstp_cond_unconverged_mask_gp->get();
      
      thrust::transform_if(
        thrust::make_constant_iterator<unsigned int>(n),
        thrust::make_constant_iterator<unsigned int>(n) + n_part,
        unconverged_mask.begin(),
        perparticle_sstp_cond.begin(),
        cuda::std::identity(),
        cuda::std::identity()
      );

      // thrust::copy_if(
      //   thrust::make_constant_iterator<unsigned int>(n),
      //   thrust::make_constant_iterator<unsigned int>(n) + n_part,
      //   unconverged_mask.begin(),
      //   perparticle_sstp_cond.begin(),
      //   cuda::std::identity()
      // );
    }
  };
};