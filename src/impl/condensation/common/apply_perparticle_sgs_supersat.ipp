namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t>
      struct add_dot_to_ssp
      { 
        const real_t dt;
        add_dot_to_ssp(const real_t dt) : dt(dt) {}

        BOOST_GPU_ENABLED
        real_t operator() (real_t ssp_val, thrust::tuple<real_t, unsigned int> dotssp_sstp_tpl) const noexcept
        {
          return ssp_val + thrust::get<0>(dotssp_sstp_tpl) * dt / thrust::get<1>(dotssp_sstp_tpl);
        }
      };
    }

    template <typename real_t, backend_t device>
    template <bool use_unconverged_mask, class it_t>
    void particles_t<real_t, device>::impl::apply_perparticle_sgs_supersat(const it_t sstp_cond_it)
    {   
      if constexpr (!use_unconverged_mask)
        thrust::transform(
          ssp.begin(), ssp.end(),  
          thrust::make_zip_iterator(thrust::make_tuple(
            dot_ssp.begin(),
            sstp_cond_it
          )),
          ssp.begin(),
          detail::add_dot_to_ssp(dt)
        );
      else
      {
        const auto &unconverged_mask = sstp_cond_unconverged_mask_gp->get();
        thrust::transform_if(
          ssp.begin(), ssp.end(),  
          thrust::make_zip_iterator(thrust::make_tuple(
            dot_ssp.begin(),
            sstp_cond_it
          )),
          unconverged_mask.begin(),
          ssp.begin(),
          detail::add_dot_to_ssp(dt),
          cuda::std::identity()  
        );
      }
    }
  };
};
