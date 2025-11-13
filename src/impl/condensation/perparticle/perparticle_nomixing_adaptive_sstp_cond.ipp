namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template<typename real_t>
      struct perparticle_nomixing_adaptive_sstp_cond_loop
      {
        const bool th_dry, const_p, turb_cond, adaptive_sstp_cond;
        const real_t dt, RH_max, cond_mlt, sstp_cond_adapt_drw2_eps, sstp_cond_adapt_drw2_max;
        const common::detail::eps_tolerance<real_t> eps_tolerance;
        const int n_dims, sstp_cond_max, sstp_cond_act;
        const RH_formula_t RH_formula;
        uintmax_t n_iter;

        perparticle_nomixing_adaptive_sstp_cond_loop(
          const opts_init_t<real_t> &opts_init,
          const opts_t<real_t> &opts,
          const int n_dims,
          const real_t &dt,
          const int &sstp_cond_max,
          const int &sstp_cond_act,
          const common::detail::eps_tolerance<real_t> &eps_tolerance, 
          const real_t &cond_mlt, 
          const uintmax_t &n_iter,
          const real_t &sstp_cond_adapt_drw2_eps,
          const real_t &sstp_cond_adapt_drw2_max
        ) : th_dry(opts_init.th_dry),
            const_p(opts_init.const_p),
            turb_cond(opts.turb_cond),
            dt(dt),
            RH_max(opts.RH_max),
            n_dims(n_dims),
            RH_formula(opts_init.RH_formula),
            eps_tolerance(eps_tolerance),
            cond_mlt(cond_mlt),
            n_iter(n_iter),
            adaptive_sstp_cond(opts_init.adaptive_sstp_cond),
            sstp_cond_act(sstp_cond_act),
            sstp_cond_max(sstp_cond_max),
            sstp_cond_adapt_drw2_eps(sstp_cond_adapt_drw2_eps),
            sstp_cond_adapt_drw2_max(sstp_cond_adapt_drw2_max)
        {}

        template<class tpl_t>
        BOOST_GPU_ENABLED void operator()(
          tpl_t tpl
        ) //noexcept
        {
          unsigned int &sstp_cond = thrust::get<0>(thrust::get<0>(tpl));
          real_t &sstp_tmp_rv = thrust::get<1>(thrust::get<0>(tpl));
          real_t &sstp_tmp_th = thrust::get<2>(thrust::get<0>(tpl));
          real_t &sstp_tmp_rh = thrust::get<3>(thrust::get<0>(tpl));
          real_t &sstp_tmp_p = thrust::get<4>(thrust::get<0>(tpl));
          const real_t &sstp_dlt_rv = thrust::get<5>(thrust::get<0>(tpl));
          const real_t &sstp_dlt_th = thrust::get<6>(thrust::get<0>(tpl));
          const real_t &sstp_dlt_rhod = thrust::get<7>(thrust::get<0>(tpl));
          const real_t &sstp_dlt_p = thrust::get<8>(thrust::get<0>(tpl));
          real_t &ssp = thrust::get<9>(thrust::get<0>(tpl));
          const real_t &dot_ssp = thrust::get<0>(thrust::get<1>(tpl));
          real_t &Tp = thrust::get<1>(thrust::get<1>(tpl));
          real_t &drw2 = thrust::get<2>(thrust::get<1>(tpl));
          real_t &drw3 = thrust::get<3>(thrust::get<1>(tpl));
          real_t &rw2 = thrust::get<4>(thrust::get<1>(tpl));
          const unsigned int &n = thrust::get<5>(thrust::get<1>(tpl));
          const real_t &dv = thrust::get<6>(thrust::get<1>(tpl));
          const real_t &lambda_D = thrust::get<7>(thrust::get<1>(tpl));
          const real_t &lambda_K = thrust::get<8>(thrust::get<1>(tpl));
          const real_t &rd3 = thrust::get<9>(thrust::get<1>(tpl));
          const real_t &kpa = thrust::get<0>(thrust::get<2>(tpl));
          const real_t &vt = thrust::get<1>(thrust::get<2>(tpl));
          const real_t &rc2 = thrust::get<2>(thrust::get<2>(tpl));

          // Two helper functions
          auto apply_noncond_perparticle_sstp_delta = [&] (const real_t &multiplier) -> void
          {
            sstp_tmp_rv += sstp_dlt_rv * multiplier;
            sstp_tmp_th += sstp_dlt_th * multiplier;
            sstp_tmp_rh += sstp_dlt_rhod * multiplier;
            if(const_p)
              sstp_tmp_p += sstp_dlt_p * multiplier;
            if(turb_cond)
              ssp += dot_ssp * dt * multiplier;
          };

          auto cond_perparticle_drw2 = [&] (
            const unsigned int &_sstp_cond,
            real_t &drw2_out
          ) -> real_t
          {
          // calculate perparticle temperature; 
            if(th_dry)
              Tp = detail::common__theta_dry__T_rhod<real_t>()(sstp_tmp_th, sstp_tmp_rh);
            else
              Tp = detail::common__theta_std__T_p<real_t>()(sstp_tmp_th, 
                thrust::make_tuple(sstp_tmp_rv, sstp_tmp_p)
              );
            // advance rw2
            if(!const_p) // sstp_tmp_p needs to be allocated even without const_p!
              sstp_tmp_p = detail::common__theta_dry__p<real_t>()(
                thrust::make_tuple(sstp_tmp_rh, sstp_tmp_rv, Tp)
              );

            const real_t RH = turb_cond ?
              detail::RH_sgs<real_t>(RH_formula)(
                thrust::make_tuple(sstp_tmp_p, sstp_tmp_rv, Tp, ssp)
              ) :
              detail::RH<real_t>(RH_formula)(
                thrust::make_tuple(sstp_tmp_p, sstp_tmp_rv, Tp)
              );
              
            drw2_out =
              detail::advance_rw2<real_t, false>(dt / _sstp_cond, RH_max, eps_tolerance, cond_mlt, n_iter)(
                rw2,
                thrust::make_tuple(
                  thrust::make_tuple(
                    sstp_tmp_rh,
                    sstp_tmp_rv,
                    Tp,
                    detail::common__vterm__visc<real_t>()(Tp),
                    rd3,
                    kpa,
                    vt,
                    lambda_D,
                    lambda_K
                  ),
                  sstp_tmp_p,
                  RH                  
                )
              ); 
            return Tp;
          };
          

          // bool converged = false;
          real_t delta_fraction_applied;
          bool first_cond_step_done_in_adaptation = sstp_cond_max == 1 ? true : false; // actually its true if sstp_cond_max is a power of 2 (?)
          // bool activates = false;

          // look for correct number of substeps
          // NOTE: this function is actually only called when adaptive_sstp_cond == true, so we skip the check below
          // if(adaptive_sstp_cond)
          {
            real_t &drw2_new = drw3; // reuse drw3 in adaptation
            real_t Tp; // temperature

            sstp_cond = sstp_cond_max; // start with max number of substeps, may be changed due to convergence or if droplets activate in this step

            // check drw convergence for increasing number of substeps
            for(int sstp_cond_try = 1; sstp_cond_try <= sstp_cond_max; sstp_cond_try*=2)
            {
              delta_fraction_applied = sstp_cond_try == 1 ? 1 : -real_t(1) / sstp_cond_try;
              apply_noncond_perparticle_sstp_delta(delta_fraction_applied);
              Tp = cond_perparticle_drw2(sstp_cond_try, sstp_cond_try == 1 ? drw2 : drw2_new);                    

              if(sstp_cond_try > 1) // check for convergence 
              {
                // static constexpr real_t tol = static_cast<real_t>(1e-9); // 1e-3 // TODO: config or opts_init parameter
                if((cuda::std::abs(drw2_new * 2 - drw2) <= sstp_cond_adapt_drw2_eps * rw2) // drw2 relative to rw2 converged
                    && cuda::std::abs(drw2 < sstp_cond_adapt_drw2_max * rw2)) // otherwise for small droplets (near activation?) drw2_new == 2*drw already for 2 substeps, but we ativate too many droplets
                // if(cuda::std::abs(drw2_new * 2 - drw2) <= tol * drw2) // drw2 converged
                {
                  // ;
                  sstp_cond = sstp_cond_try / 2;
                  apply_noncond_perparticle_sstp_delta(-delta_fraction_applied); // revert last addition to get to a state after one step of converged number            
                  first_cond_step_done_in_adaptation = true;
                  // first_cond_step_done_in_adaptation = false;
                  // fprintf(stderr, "rd3: %g converged with sstp_cond = %u drw2: %g drw2_new: %g rw2: %g\n", rd3, sstp_cond, drw2, drw2_new, rw2);
                  break;
                }
                drw2 = drw2_new;                
              }
            }

            // override number of substeps for SDs that de/activate in this timestep;
            if(sstp_cond_act > 1)
            // --> pimpl->set_activating_perparticle_sstp_cond(pimpl->sstp_cond_act); 
            {
              // real_t &rc2 = drw3; // reuse drw3 as it is still not needed here
              // TODO: how much does rc depend on T? maybe skip this dependance and calc rc2 only if rd changes?
              // rc2 = detail::rw3_cr<real_t>()(rd3, thrust::make_tuple(kpa, Tp));
              // rc2 = 1;


              if ( ( rw2 < rc2 && (rw2 + sstp_cond * drw2) > rc2 ) || 
                   ( rw2 > rc2 && (rw2 + sstp_cond * drw2) < rc2 ) )
              {
                sstp_cond = sstp_cond_act;
                first_cond_step_done_in_adaptation = false;
                // apply_noncond_perparticle_sstp_delta(delta_fraction_applied); // revert to state before adaptation loop
                // activates = true; 
              }
            }
            // TODO: sstp_tmp and drw2 are not correct for these!

            // if(first_cond_step_done_in_adaptation)
            // {
            //   // apply_noncond_perparticle_sstp_delta(delta_fraction_applied); // revert to state before adaptation loop (beacause sstp_cond == sstp_cond_max and sstp_cond_max may not be a power of 2)
            //   apply_noncond_perparticle_sstp_delta(-delta_fraction_applied); // revert last addition to get to a state after one step of converged number            
            // }
            // else
            if(!first_cond_step_done_in_adaptation)
            {
              apply_noncond_perparticle_sstp_delta(delta_fraction_applied); // revert to state before adaptation loop (beacause sstp_cond == sstp_cond_max and sstp_cond_max may not be a power of 2)
            }

            // pimpl->reset_perparticle_sstp_tmp_and_ssp_before_substepping(); // TODO: reuse it and not reset?
          
//           // in drw2, we have correct drw2 for each particle; in sstp_count we have number of substeps needed;
//           // in the first subsequent step we can skip cond_perparticle_rw2_change for converged particles, and add stored drw3 (or make it drw2?) instead
//           // sstp_tmp and ssp are one step too far ?

          }            

          delta_fraction_applied = real_t(1) / sstp_cond;
          // actual condensation substepping
          for(int step = 0; step < sstp_cond; ++step)
          {
            // ---> pimpl->apply_noncond_perparticle_sstp_delta(real_t(1) / pimpl->sstp_cond);
            // This idea wont work if sstp_cond is not a power of 2!!! then if notconverged, sstp_cond can be different than sstp_cond_max (?)
            if(!first_cond_step_done_in_adaptation || step > 0) // || !adaptive_sstp_cond) // latter part would be required if we used this function also when not adapting (possible future extension)
              apply_noncond_perparticle_sstp_delta(delta_fraction_applied);
            // sstp_tmp_rv += sstp_dlt_rv * (real_t(1) / sstp_cond);
            // sstp_tmp_th += sstp_dlt_th * (real_t(1) / sstp_cond);
            // sstp_tmp_rh += sstp_dlt_rhod * (real_t(1) / sstp_cond);
            // if(const_p)
            //   sstp_tmp_p += sstp_dlt_p * (real_t(1) / sstp_cond);
              
            // if(turb_cond)
            // ---> pimpl->apply_perparticle_sgs_supersat(pimpl->dt / pimpl->sstp_cond);
              // ssp += dot_ssp * dt * (real_t(1) / sstp_cond);

            // ---> pimpl->cond_perparticle_drw2(pimpl->dt / pimpl->sstp_cond, opts.RH_max, opts.turb_cond, pimpl->drw2_gp->get());
            if(!first_cond_step_done_in_adaptation || step > 0) // || !adaptive_sstp_cond) // same as above
              cond_perparticle_drw2(sstp_cond, drw2);
            
            // ---> pimpl->cond_perparticle_drw3_from_drw2();
            drw3 = detail::drw2_to_drw3<real_t>()(drw2, rw2);

            // ---> pimpl->apply_perparticle_drw2();
            rw2 += drw2;

            // ---> pimpl->apply_perparticle_drw3_to_perparticle_rv_and_th();
            drw3 = detail::rw3diff2drv<real_t>(
              - common::moist_air::rho_w<real_t>() / si::kilograms * si::cubic_metres
              * real_t(4./3) * real_t(3.14159265358979323846264338), n_dims
            ) (drw3, 
              thrust::make_tuple(sstp_tmp_rh, n, dv)
            );

            sstp_tmp_rv += drw3;

            drw3 = detail::dth<real_t>()(
              thrust::make_tuple(drw3, Tp, sstp_tmp_th)
            );

            sstp_tmp_th += drw3;
          }
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::perparticle_nomixing_adaptive_sstp_cond(const opts_t<real_t> &opts) { 

      auto &perparticle_sstp_cond = perparticle_sstp_cond_gp->get();
      auto &sstp_dlt_rv = sstp_dlt_rv_gp->get();
      auto &sstp_dlt_th = sstp_dlt_th_gp->get();
      auto &sstp_dlt_rhod = sstp_dlt_rhod_gp->get();
      auto &sstp_dlt_p = sstp_dlt_p_gp->get();
      auto &Tp = Tp_gp->get();
      auto &drw2 = drw2_gp->get();
      auto &drw3 = drw3_gp->get();
      const auto &lambda_D = lambda_D_gp->get(); 
      const auto &lambda_K = lambda_K_gp->get(); 
      
      auto pptcl_nomix_sstp_cond_args_zip = 
        thrust::make_zip_iterator(thrust::make_tuple(
          thrust::make_zip_iterator(thrust::make_tuple(
            perparticle_sstp_cond.begin(),
            sstp_tmp_rv.begin(),
            sstp_tmp_th.begin(),
            sstp_tmp_rh.begin(),
            sstp_tmp_p.begin(),
            sstp_dlt_rv.begin(),
            sstp_dlt_th.begin(),
            sstp_dlt_rhod.begin(),
            sstp_dlt_p.begin(),
            ssp.begin()
          )),
          thrust::make_zip_iterator(thrust::make_tuple(
            dot_ssp.begin(),
            Tp.begin(),
            drw2.begin(),
            drw3.begin(),
            rw2.begin(),
            n.begin(),
            thrust::make_permutation_iterator(dv.begin(), ijk.begin()),
            thrust::make_permutation_iterator(lambda_D.begin(), ijk.begin()),
            thrust::make_permutation_iterator(lambda_K.begin(), ijk.begin()),
            rd3.begin()
          )),
          thrust::make_zip_iterator(thrust::make_tuple(
            kpa.begin(),
            vt.begin(),
            rc2.begin()
          ))
        ));

      thrust::for_each(
        pptcl_nomix_sstp_cond_args_zip,
        pptcl_nomix_sstp_cond_args_zip + n_part,
        detail::perparticle_nomixing_adaptive_sstp_cond_loop<real_t>(
          opts_init, opts, n_dims, dt, sstp_cond, sstp_cond_act, config.eps_tolerance, config.cond_mlt, config.n_iter, config.sstp_cond_adapt_drw2_eps, config.sstp_cond_adapt_drw2_max
        )
      );
    }
  };
};