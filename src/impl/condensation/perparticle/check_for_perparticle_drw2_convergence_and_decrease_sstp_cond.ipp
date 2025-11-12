namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template<typename real_t>
      struct convergence_checker
      {
        // static constexpr real_t tol = static_cast<real_t>(1e-6); // 1e-3 // TODO: config or opts_init parameter
        real_t dt_ratio, sstp_cond_adapt_eps;
        convergence_checker(const real_t dt_ratio_, const real_t sstp_cond_adapt_eps_) : dt_ratio(dt_ratio_), sstp_cond_adapt_eps(sstp_cond_adapt_eps_) {}

        template<class tpl_t>
        BOOST_GPU_ENABLED bool operator()(unsigned int &sstp_cond, tpl_t rw2_drw2_new_old) noexcept
        {
          const auto &rw2 = cuda::std::get<0>(rw2_drw2_new_old);
          const auto &drw2_new = cuda::std::get<1>(rw2_drw2_new_old); // drw2 for new (larger) number of substeps
          const auto &drw2_old = cuda::std::get<2>(rw2_drw2_new_old); // drw2 for old (smaller) number of substeps
          if
          (
            // cuda::std::abs(drw2_new * dt_ratio - drw2_old) > tol * cuda::std::abs(drw2_old) || // drw2 not converged
            cuda::std::abs(drw2_new * dt_ratio - drw2_old) > sstp_cond_adapt_eps * rw2 // drw2 relative to rw2 not converged
            // TODO: check for drw2 < const * rw2, as in perparticle_nomixing_adaptive_sstp_cond_loop
          ) 
            return true; // not converged
          else // converged
          {
            if(sstp_cond > 1)
              sstp_cond /= dt_ratio;
            return false;
          }
        }
      };

    };
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::check_for_perparticle_drw2_convergence_and_decrease_sstp_cond(
      const thrust_device::vector<real_t> &drw2,
      thrust_device::vector<real_t> &drw2_old,
      const real_t dt_ratio
    ) { 
      auto &unconverged_mask = sstp_cond_unconverged_mask_gp->get();
      auto &perparticle_sstp_cond = perparticle_sstp_cond_gp->get();
      thrust::transform_if(
        perparticle_sstp_cond.begin(), perparticle_sstp_cond.end(),
        thrust::make_zip_iterator(thrust::make_tuple(
          rw2.begin(),
          drw2.begin(),
          drw2_old.begin()
        )),
        unconverged_mask.begin(),
        unconverged_mask.begin(),
        detail::convergence_checker<real_t>(dt_ratio, config.sstp_cond_adapt_eps),
        cuda::std::identity()
      );
    }
  };
};