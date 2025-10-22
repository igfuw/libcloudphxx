
namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::apply_perparticle_sgs_supersat(
    const real_t &dt
    )
    {   
    namespace arg = thrust::placeholders;

    thrust::transform(
        ssp.begin(), ssp.end(),  
        dot_ssp.begin(),
        ssp.begin(),
        arg::_1 + arg::_2 * dt
    );
    }
  };
};