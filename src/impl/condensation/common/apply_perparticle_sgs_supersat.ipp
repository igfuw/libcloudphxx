namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    template <bool use_unconverged_mask>
    void particles_t<real_t, device>::impl::apply_perparticle_sgs_supersat(const real_t &dt) 
    {   
      namespace arg = thrust::placeholders;

      if constexpr (!use_unconverged_mask)
        thrust::transform(
          ssp.begin(), ssp.end(),  
          dot_ssp.begin(),
          ssp.begin(),
          arg::_1 + arg::_2 * dt
        );
      else
      {
        const auto &unconverged_mask = cond_sstp_unconverged_mask_gp->get();
        thrust::transform_if(
            ssp.begin(), ssp.end(),  
            dot_ssp.begin(),
            unconverged_mask.begin(),
            ssp.begin(),
            arg::_1 + arg::_2 * dt,
            cuda::std::identity()  
        );
      }
    }
  };
};
