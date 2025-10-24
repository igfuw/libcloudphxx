namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::set_perparticle_cond_sstp(const unsigned int &n) noexcept
    {
      auto &perparticle_cond_sstp = perparticle_cond_sstp_gp->get();

      thrust::fill(
        perparticle_cond_sstp.begin(),
        perparticle_cond_sstp.end(),
        n
      );
    }
  };
};