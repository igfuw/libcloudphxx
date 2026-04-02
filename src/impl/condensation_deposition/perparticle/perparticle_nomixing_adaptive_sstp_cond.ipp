namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template<typename real_t>
      struct perparticle_nomixing_adaptive_sstp_cond_loop
      {
        const bool th_dry, const_p, turb_cond, adaptive_sstp_cond, ice_switch, cond, depo;
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
          const uintmax_t &n_iter
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
            sstp_cond_adapt_drw2_eps(opts_init.sstp_cond_adapt_drw2_eps),
            sstp_cond_adapt_drw2_max(opts_init.sstp_cond_adapt_drw2_max),
            ice_switch(opts_init.ice_switch),
            cond(opts.cond),
            depo(opts.depo)
        {}

        template<class tpl_t>
        BOOST_GPU_ENABLED void operator()(
          tpl_t tpl
        ) //noexcept
        {
          // copy values into local variables
          // variables that are not modified
          const real_t sstp_dlt_rv = thrust::get<5>(thrust::get<0>(tpl));
          const real_t sstp_dlt_th = thrust::get<6>(thrust::get<0>(tpl));
          const real_t sstp_dlt_rhod = thrust::get<7>(thrust::get<0>(tpl));
          const real_t sstp_dlt_p = thrust::get<8>(thrust::get<0>(tpl));
          const auto n = thrust::get<2>(thrust::get<1>(tpl));
          const real_t dv = thrust::get<3>(thrust::get<1>(tpl));
          const real_t lambda_D = thrust::get<4>(thrust::get<1>(tpl));
          const real_t lambda_K = thrust::get<5>(thrust::get<1>(tpl));
          const real_t rd3 = thrust::get<6>(thrust::get<1>(tpl));
          const real_t kpa = thrust::get<0>(thrust::get<2>(tpl));
          const real_t vt = thrust::get<1>(thrust::get<2>(tpl));
          const real_t dot_ssp = turb_cond ? thrust::get<0>(thrust::get<1>(tpl)) : 0;

          // variables that are modified, we make local copies regardless and copy back at the end
          unsigned int sstp_cond; // its set in this function, old value not important
          real_t sstp_tmp_rv = thrust::get<1>(thrust::get<0>(tpl));
          real_t sstp_tmp_th = thrust::get<2>(thrust::get<0>(tpl));
          real_t sstp_tmp_rh = thrust::get<3>(thrust::get<0>(tpl));
          real_t sstp_tmp_p = const_p ? thrust::get<4>(thrust::get<0>(tpl)) : 0;
          real_t ssp = turb_cond ? thrust::get<9>(thrust::get<0>(tpl)) : 0;
          real_t rw2 = thrust::get<1>(thrust::get<1>(tpl));
          thrust::tuple<real_t, real_t> ice_ac; // = thrust::make_tuple(0, 0); // only used if ice_switch is true, but we need to define it here to be able to use it in lambdas
          real_t ice_rho; 

          if(rw2 > 0 && !cond) return; // skip liquid droplets if condensation is turned off

          if(ice_switch)
          {
            real_t ice_a = thrust::get<7>(thrust::get<1>(tpl));
            if(ice_a > 0 && !depo) return; // skip ice particles if deposition is turned off
            real_t ice_c = thrust::get<8>(thrust::get<1>(tpl));
            ice_ac = thrust::make_tuple(ice_a, ice_c);
            ice_rho = thrust::get<9>(thrust::get<1>(tpl));
          }

          const bool ice = rw2<=0; // flag to indicate if the particle is ice or liquid;

          real_t drw2, Tp, RH;
          thrust::tuple<real_t, real_t> dice_ac;

          // helper functions
          auto _apply_noncond_perparticle_sstp_delta = [&] (const real_t &multiplier) -> void
          {
            sstp_tmp_rv += sstp_dlt_rv * multiplier;
            sstp_tmp_th += sstp_dlt_th * multiplier;
            sstp_tmp_rh += sstp_dlt_rhod * multiplier;
            if(const_p)
              sstp_tmp_p += sstp_dlt_p * multiplier;
            if(turb_cond)
              ssp += dot_ssp * dt * multiplier;
          };

          auto _calc_Tp = [&] () -> void
          {
            if(th_dry)
              Tp = detail::common__theta_dry__T_rhod<real_t>()(sstp_tmp_th, sstp_tmp_rh);
            else
              Tp = detail::common__theta_std__T_p<real_t>()(sstp_tmp_th, 
                thrust::make_tuple(sstp_tmp_rv, sstp_tmp_p)
              );
          };

          auto _calc_sstp_tmp_p = [&] () -> void
          {
            if(!const_p) // sstp_tmp_p needs to be allocated even without const_p!
              sstp_tmp_p = detail::common__theta_dry__p<real_t>()(
                thrust::make_tuple(sstp_tmp_rh, sstp_tmp_rv, Tp)
              );
          };

          auto _calc_RH = [&] () -> void
          {
            if(!ice)
              RH = turb_cond ?
                detail::RH_hlpr<real_t, false, true>(RH_formula)(
                  thrust::make_tuple(sstp_tmp_p, sstp_tmp_rv, Tp, ssp)
                ) :
                detail::RH_hlpr<real_t, false, false>(RH_formula)(
                  thrust::make_tuple(sstp_tmp_p, sstp_tmp_rv, Tp)
                );
            else
              RH = turb_cond ?
                detail::RH_hlpr<real_t, true, true>(RH_formula)(
                  thrust::make_tuple(sstp_tmp_p, sstp_tmp_rv, Tp, ssp)
                ) :
                detail::RH_hlpr<real_t, true, false>(RH_formula)(
                  thrust::make_tuple(sstp_tmp_p, sstp_tmp_rv, Tp)
                );
          };

          // bool converged = false;
          real_t delta_fraction_applied;
          bool first_cond_step_done_in_adaptation = sstp_cond_max == 1 ? true : false; // actually its true if sstp_cond_max is a power of 2 (?)
          // bool activates = false;

          // look for correct number of substeps
          // NOTE: this function is actually only called when adaptive_sstp_cond == true, so we skip the check below
          // if(adaptive_sstp_cond)
          {
            real_t drw2_new; 
            thrust::tuple<real_t, real_t> dice_ac_new;
            // real_t Tp; // temperature

            sstp_cond = sstp_cond_max; // start with max number of substeps, may be changed due to convergence or if droplets activate in this step

            // check dX convergence for increasing number of substeps
            for(int sstp_cond_try = 1; sstp_cond_try <= sstp_cond_max; sstp_cond_try*=2)
            {
              delta_fraction_applied = sstp_cond_try == 1 ? 1 : -real_t(1) / sstp_cond_try;
              _apply_noncond_perparticle_sstp_delta(delta_fraction_applied);
              // _cond_perparticle_drw2(sstp_cond_try, sstp_cond_try == 1 ? drw2 : drw2_new); // also updates Tp!

              _calc_Tp();              
              _calc_sstp_tmp_p();
              _calc_RH();
              if(!ice)
                (sstp_cond_try == 1 ? drw2 : drw2_new) = 
                    detail::advance_rw2<real_t, false>(dt / sstp_cond_try, RH_max, eps_tolerance, cond_mlt, n_iter)(
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
              else // ice
                (sstp_cond_try == 1 ? dice_ac : dice_ac_new) = 
                  detail::advance_ice_ac<real_t, false>(dt / sstp_cond_try, RH_max)(
                    ice_ac,
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

              if(sstp_cond_try > 1) // check for convergence 
              {
                // TODO: get rid of _max and of actviation adaptation?
                if(!ice)
                {
                  if((cuda::std::abs(drw2_new * 2 - drw2) <= sstp_cond_adapt_drw2_eps * rw2) // drw2 relative to rw2 converged
                      && cuda::std::abs(drw2 < sstp_cond_adapt_drw2_max * rw2)) // otherwise for small droplets (near activation?) drw2_new == 2*drw already for 2 substeps, but we ativate too many droplets
                  // if(cuda::std::abs(drw2_new * 2 - drw2) <= tol * drw2) // drw2 converged
                  {
                    sstp_cond = sstp_cond_try / 2;
                    _apply_noncond_perparticle_sstp_delta(-delta_fraction_applied); // revert last addition to get to a state after one step of converged number            
                    first_cond_step_done_in_adaptation = true;
                    break;
                  }
                  drw2 = drw2_new;                
                }
                else // ice
                {
                  const real_t da_dt_rel = cuda::std::abs(thrust::get<0>(dice_ac_new) * 2 - thrust::get<0>(dice_ac)) / (thrust::get<0>(dice_ac) + 1e-20);
                  const real_t dc_dt_rel = cuda::std::abs(thrust::get<1>(dice_ac_new) * 2 - thrust::get<1>(dice_ac)) / (thrust::get<1>(dice_ac) + 1e-20);
                  if(da_dt_rel <= sstp_cond_adapt_drw2_eps && da_dt_rel <= sstp_cond_adapt_drw2_eps
                      && da_dt_rel <= sstp_cond_adapt_drw2_max && dc_dt_rel <= sstp_cond_adapt_drw2_max)
                  {
                    sstp_cond = sstp_cond_try / 2;
                    _apply_noncond_perparticle_sstp_delta(-delta_fraction_applied); // revert last addition to get to a state after one step of converged number            
                    first_cond_step_done_in_adaptation = true;
                    break;
                  }
                  dice_ac = dice_ac_new;
                }
              }
            }

            // override number of substeps for SDs that de/activate in this timestep;
            if(sstp_cond_act > 1)
            {
              const real_t rc2 = thrust::get<2>(thrust::get<2>(tpl));

              if ( ( rw2 < rc2 && (rw2 + sstp_cond * drw2) > rc2 ) || 
                   ( rw2 > rc2 && (rw2 + sstp_cond * drw2) < rc2 ) )
              {
                sstp_cond = sstp_cond_act;
                first_cond_step_done_in_adaptation = false;
              }
            }
            if(!first_cond_step_done_in_adaptation)
            {
              _apply_noncond_perparticle_sstp_delta(sstp_cond_max == 1 ? -delta_fraction_applied : delta_fraction_applied); // revert to state before adaptation loop (beacause sstp_cond == sstp_cond_max and sstp_cond_max may not be a power of 2); If only one step was tried, whole change was applied; If more steps were tried, we are moving back from the entire step
            }
          }            

          delta_fraction_applied = real_t(1) / sstp_cond;

          // actual condensation substepping
          if(!ice)
          {
            auto _advance_rw2 = detail::advance_rw2<real_t, true>(dt / sstp_cond, RH_max, eps_tolerance, cond_mlt, n_iter);
            real_t &rw3 = drw2; // drw2 needed only at the start of the first step
            real_t drw3;

            auto rw2torw3 = detail::rw2torwX<real_t, 3>();

            for(int step = 0; step < sstp_cond; ++step)
            {
              drw3 = step > 0 ? -rw3 : -rw2torw3(rw2);

              if(first_cond_step_done_in_adaptation && step == 0)
              {
                rw2 += drw2;
              }
              else
              {
                _apply_noncond_perparticle_sstp_delta(delta_fraction_applied);
                _calc_Tp();              
                _calc_sstp_tmp_p();
                _calc_RH();

                rw2 = _advance_rw2(
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
              }

              if (step < sstp_cond - 1)
              {
                rw3 = rw2torw3(rw2);
                drw3 += rw3;
              }
              else
                drw3 += rw2torw3(rw2);
              
              drw3 = detail::massdiff2drv<real_t>(
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
          else // actual deposition substepping
          {
            auto _advance_ice_ac = detail::advance_ice_ac<real_t, true>(dt / sstp_cond, RH_max);
            real_t &ice_mass = drw2; // if its ice, drw2 is not used at all
            real_t dice_mass; 

            auto ice_mass_fctr = detail::ice_mass<real_t>();

            for(int step = 0; step < sstp_cond; ++step)
            {
              dice_mass = step > 0 ? -ice_mass : -ice_mass_fctr(thrust::make_tuple(thrust::get<0>(ice_ac), thrust::get<1>(ice_ac), ice_rho));

              if(first_cond_step_done_in_adaptation && step == 0)
              {
                thrust::get<0>(ice_ac) += thrust::get<0>(dice_ac);
                thrust::get<1>(ice_ac) += thrust::get<1>(dice_ac);
              }
              else
              {
                _apply_noncond_perparticle_sstp_delta(delta_fraction_applied);
                _calc_Tp();              
                _calc_sstp_tmp_p();
                _calc_RH();

                ice_ac = _advance_ice_ac(
                  ice_ac,
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
              }

              if (step < sstp_cond - 1)
              {
                ice_mass = ice_mass_fctr(thrust::make_tuple(thrust::get<0>(ice_ac), thrust::get<1>(ice_ac), ice_rho));
                dice_mass += ice_mass;
              }
              else
                dice_mass += ice_mass_fctr(thrust::make_tuple(thrust::get<0>(ice_ac), thrust::get<1>(ice_ac), ice_rho));
              
              dice_mass = detail::massdiff2drv<real_t>(real_t(1), n_dims) 
                (dice_mass, 
                  thrust::make_tuple(sstp_tmp_rh, n, dv)
                );

              sstp_tmp_rv += dice_mass;

              dice_mass = detail::dth_dep<real_t>()(
                thrust::make_tuple(dice_mass, Tp, sstp_tmp_th)
              );

              sstp_tmp_th += dice_mass;
            }
          }

          // copy back modified variables
          thrust::get<0>(thrust::get<0>(tpl)) = sstp_cond;
          thrust::get<1>(thrust::get<0>(tpl)) = sstp_tmp_rv;
          thrust::get<2>(thrust::get<0>(tpl)) = sstp_tmp_th;
          thrust::get<3>(thrust::get<0>(tpl)) = sstp_tmp_rh;
          if(const_p)
            thrust::get<4>(thrust::get<0>(tpl)) = sstp_tmp_p;
          if(turb_cond)
            thrust::get<9>(thrust::get<0>(tpl)) = ssp;   
          if(!ice)
            thrust::get<1>(thrust::get<1>(tpl)) = rw2;
          else // ice
          {
            thrust::get<7>(thrust::get<1>(tpl)) = thrust::get<0>(ice_ac);
            thrust::get<8>(thrust::get<1>(tpl)) = thrust::get<1>(ice_ac);
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
      // auto &Tp = Tp_gp->get();
      // auto &drwX = drwX_gp->get();
      // auto &rwX = rwX_gp->get();
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
            // Tp.begin(),
            // drwX.begin(),
            // rwX.begin(),
            rw2.begin(),
            n.begin(),
            thrust::make_permutation_iterator(dv.begin(), ijk.begin()),
            thrust::make_permutation_iterator(lambda_D.begin(), ijk.begin()),
            thrust::make_permutation_iterator(lambda_K.begin(), ijk.begin()),
            rd3.begin(),
            ice_a.begin(),
            ice_c.begin(),
            ice_rho.begin()
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
          opts_init, opts, n_dims, dt, sstp_cond, sstp_cond_act, config.eps_tolerance, config.cond_mlt, config.n_iter
        )
      );
    }
  };
};
