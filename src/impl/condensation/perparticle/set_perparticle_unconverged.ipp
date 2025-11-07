namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::set_perparticle_unconverged() noexcept
    {
      auto &unconverged_mask = sstp_cond_unconverged_mask_gp->get();
      thrust::fill(
        unconverged_mask.begin(),
        unconverged_mask.end(),
        true
      );
    }
  };
};