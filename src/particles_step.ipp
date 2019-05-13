// vim:filetype=cpp

/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief timestepping routine for super droplets
  */

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::step_sync(
      const opts_t<real_t> &opts,
      arrinfo_t<real_t> th,
      arrinfo_t<real_t> rv,
      const arrinfo_t<real_t> rhod,      // defaults to NULL-NULL pair (e.g. kinematic or boussinesq model)
      const arrinfo_t<real_t> courant_x, // defaults to NULL-NULL pair (e.g. kinematic model)
      const arrinfo_t<real_t> courant_y, // defaults to NULL-NULL pair (e.g. kinematic model)
      const arrinfo_t<real_t> courant_z, // defaults to NULL-NULL pair (e.g. kinematic model)
      std::map<enum chem_species_t, arrinfo_t<real_t> > ambient_chem
    )
    {
      sync_in(th, rv, rhod, courant_x, courant_y, courant_z, ambient_chem);
      step_cond(opts, th, rv, ambient_chem);
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::sync_in(
      arrinfo_t<real_t> th,
      arrinfo_t<real_t> rv,
      const arrinfo_t<real_t> rhod,      // defaults to NULL-NULL pair (e.g. kinematic or boussinesq model)
      const arrinfo_t<real_t> courant_x, // defaults to NULL-NULL pair (e.g. kinematic model)
      const arrinfo_t<real_t> courant_y, // defaults to NULL-NULL pair (e.g. kinematic model)
      const arrinfo_t<real_t> courant_z, // defaults to NULL-NULL pair (e.g. kinematic model)
      std::map<enum chem_species_t, arrinfo_t<real_t> > ambient_chem
    )
    {
      // sanity checks
      if (!pimpl->init_called)
        throw std::runtime_error("please call init() before calling step_sync()");
      if (pimpl->should_now_run_async)
        throw std::runtime_error("please call step_async() before calling step_sync() again");

      if (th.is_null() || rv.is_null())
        throw std::runtime_error("passing th and rv is mandatory");

 // <TODO> - code duplicated from init() !
      if (!courant_x.is_null() || !courant_y.is_null() || !courant_z.is_null())
      {
	if (pimpl->n_dims == 0)
	  throw std::runtime_error("Courant numbers passed in 0D setup");

	if (pimpl->n_dims == 1 && (courant_x.is_null() || !courant_y.is_null() || !courant_z.is_null()))
	  throw std::runtime_error("Only X Courant number allowed in 1D setup");

	if (pimpl->n_dims == 2 && (courant_x.is_null() || !courant_y.is_null() || courant_z.is_null()))
	  throw std::runtime_error("Only X and Z Courant numbers allowed in 2D setup");

	if (pimpl->n_dims == 3 && (courant_x.is_null() || courant_y.is_null() || courant_z.is_null()))
	  throw std::runtime_error("All XYZ Courant number components required in 3D setup");
      }

      if (pimpl->opts_init.chem_switch && ambient_chem.size() != chem_gas_n)
        throw std::runtime_error("chemistry was not switched off and ambient_chem is empty");

      if (!pimpl->opts_init.chem_switch && ambient_chem.size() != 0)
        throw std::runtime_error("chemistry was switched off and ambient_chem is not empty");
// </TODO>

      if (pimpl->l2e[&pimpl->courant_x].size() == 0) // TODO: y, z,...
      {
        // TODO: many max or m1 used, unify it
#if !defined(__NVCC__)
        using std::max;
#endif
        // TODO: copy-pasted from init
        if (!courant_x.is_null()) pimpl->init_e2l(courant_x, &pimpl->courant_x, 1, 0, 0, - pimpl->halo_x);
        if (!courant_y.is_null()) pimpl->init_e2l(courant_y, &pimpl->courant_y, 0, 1, 0, pimpl->n_x_bfr * pimpl->opts_init.nz - pimpl->halo_y);
        if (!courant_z.is_null()) pimpl->init_e2l(courant_z, &pimpl->courant_z, 0, 0, 1, pimpl->n_x_bfr * max(1, pimpl->opts_init.ny) - pimpl->halo_z);
      }

      // did rhod change
      pimpl->var_rho = !rhod.is_null();

      // syncing in Eulerian fields (if not null)
      pimpl->sync(th,             pimpl->th);
      pimpl->sync(rv,             pimpl->rv);
      pimpl->sync(courant_x,      pimpl->courant_x);
      pimpl->sync(courant_y,      pimpl->courant_y);
      pimpl->sync(courant_z,      pimpl->courant_z);
      pimpl->sync(rhod,           pimpl->rhod);

      nancheck(pimpl->th, " th after sync-in");
      nancheck(pimpl->rv, " rv after sync-in");
      nancheck(pimpl->courant_x, " courant_x after sync-in");
      nancheck(pimpl->courant_y, " courant_y after sync-in");
      nancheck(pimpl->courant_z, " courant_z after sync-in");
      nancheck(pimpl->rhod, " rhod after sync-in");

      assert(*thrust::min_element(pimpl->rv.begin(), pimpl->rv.end()) >= 0);
      assert(*thrust::min_element(pimpl->th.begin(), pimpl->th.end()) >= 0);
      assert(*thrust::min_element(pimpl->rhod.begin(), pimpl->rhod.end()) >= 0);

      // check if courants are greater than 2 since it would break the predictor-corrector (halo of size 2 in the x direction) 
        if(!(pimpl->opts_init.adve_scheme != as_t::pred_corr || (courant_x.is_null() || ((*(thrust::min_element(pimpl->courant_x.begin(), pimpl->courant_x.end()))) >= real_t(-2.) )) ))
        {
#if !defined(NDEBUG)
          std::cerr << "Courant x less than -2 in pred_corr advection, minimum courant x: " << *(thrust::min_element(pimpl->courant_x.begin(), pimpl->courant_x.end())) << " at " << (thrust::min_element(pimpl->courant_x.begin(), pimpl->courant_x.end())) - pimpl->courant_x.begin() << " using first order advection in this step" << std::endl;
#endif
          pimpl->adve_scheme = as_t::euler;
        }
        if(!(pimpl->opts_init.adve_scheme != as_t::pred_corr || (courant_x.is_null() || ((*(thrust::max_element(pimpl->courant_x.begin(), pimpl->courant_x.end()))) <= real_t(2.) )) ))
        {
#if !defined(NDEBUG)
          std::cerr << "Courant x more than 2 in pred_corr advection, maximum courant x: " << *(thrust::max_element(pimpl->courant_x.begin(), pimpl->courant_x.end())) << " at " << (thrust::max_element(pimpl->courant_x.begin(), pimpl->courant_x.end())) - pimpl->courant_x.begin() <<  " using first order advection in this step" << std::endl;
#endif
          pimpl->adve_scheme = as_t::euler;
        }

      if (pimpl->opts_init.chem_switch){
        for (int i = 0; i < chem_gas_n; ++i){
          pimpl->sync(
            ambient_chem.at((chem_species_t)i), 
            pimpl->ambient_chem[(chem_species_t)i]
          );
        }
      }
  
      pimpl->should_now_run_cond = true;
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::step_cond(
      const opts_t<real_t> &opts,
      arrinfo_t<real_t> th,                                            // for sync-out
      arrinfo_t<real_t> rv,                                            // for sync-out
      std::map<enum chem_species_t, arrinfo_t<real_t> > ambient_chem   // for sync-out
    )
    {
      //sanity checks
      if (!pimpl->should_now_run_cond)
        throw std::runtime_error("please call sync_in() before calling step_cond()");

      pimpl->should_now_run_cond = false;

      // condensation/evaporation 
      if (opts.cond) 
      {
        if(pimpl->opts_init.exact_sstp_cond && pimpl->opts_init.sstp_cond > 1)
        // apply substeps per-particle logic
        {
          for (int step = 0; step < pimpl->opts_init.sstp_cond; ++step) 
          {   
            pimpl->sstp_step_exact(step);
            pimpl->cond_sstp(pimpl->opts_init.dt / pimpl->opts_init.sstp_cond, opts.RH_max); 
          } 
          // copy sstp_tmp_rv and th to rv and th
          pimpl->update_state(pimpl->rv, pimpl->sstp_tmp_rv);
          pimpl->update_state(pimpl->th, pimpl->sstp_tmp_th);
        }
        else
        // apply per-cell sstp logic
        {
          for (int step = 0; step < pimpl->opts_init.sstp_cond; ++step) 
          {   
            pimpl->sstp_step(step);
            pimpl->hskpng_Tpr(); 
            pimpl->cond(pimpl->opts_init.dt / pimpl->opts_init.sstp_cond, opts.RH_max);
          }
        }

        nancheck(pimpl->th, " th after cond");
        nancheck(pimpl->rv, " rv after cond");

        // saving rv to be used as rv_old
        pimpl->sstp_save();
      }

      // chemistry
      // TODO: chemistry substepping still done the old way, i.e. per cell not per particle
      // TODO2: shouldn't we run hskpng_Tpr before chemistry?
      if (opts.chem_dsl or opts.chem_dsc or opts.chem_rct) 
      {
        for (int step = 0; step < pimpl->opts_init.sstp_chem; ++step) 
        {   
          // calculate new volume of droplets (needed for chemistry)
          pimpl->chem_vol_ante();

          // set flag for those SD that are big enough to have chemical reactions
          pimpl->chem_flag_ante();

          if (opts.chem_dsl)
          {
            //adjust trace gases to substepping
            pimpl->sstp_step_chem(step);

            //dissolving trace gases (Henrys law)
            pimpl->chem_henry(pimpl->opts_init.dt / pimpl->opts_init.sstp_chem);

            //cleanup - TODO think of something better
            pimpl->chem_cleanup();
          }
 
          if (opts.chem_dsc)
          { //dissociation
            pimpl->chem_dissoc();

            //cleanup - TODO think of something better
            pimpl->chem_cleanup();
          }
           
          if (opts.chem_rct)
          { //oxidation 
            pimpl->chem_react(pimpl->opts_init.dt / pimpl->opts_init.sstp_chem);

            //cleanup - TODO think of something better
            pimpl->chem_cleanup();
          }
        }
      }

      // aerosol source, in sync since it changes th/rv
      if (opts.src) 
      {
        // sanity check
        if (pimpl->opts_init.src_switch == false) throw std::runtime_error("aerosol source was switched off in opts_init");

        // update the step counter since src was turned on
        ++pimpl->stp_ctr;

        // introduce new particles with the given time interval
        if(pimpl->stp_ctr == pimpl->opts_init.supstp_src) 
        {
          pimpl->src(pimpl->opts_init.supstp_src * pimpl->opts_init.dt);
        }
      }
      else pimpl->stp_ctr = 0; //reset the counter if source was turned off

      if(opts.cond || pimpl->stp_ctr == pimpl->opts_init.supstp_src)
      {
        // syncing out // TODO: this is not necesarry in off-line mode (see coupling with DALES)
        pimpl->sync(pimpl->th, th);
        pimpl->sync(pimpl->rv, rv);
        pimpl->stp_ctr = 0; //reset the counter
      }

      if (opts.chem_dsl == true)
      {
        // syncing out trace gases // TODO: this is not necesarry in off-line mode (see coupling with DALES)
        for (int i = 0; i < chem_gas_n; ++i)
          pimpl->sync(
            pimpl->ambient_chem[(chem_species_t)i],
            ambient_chem.at((chem_species_t)i)
          );
      }

      // TODO: revert to old th and rv, so that processes in step_async see the same th and rv as passed in sync_in?;
      //       but then what happens in step_cond after condensation (i.e. chemistry and source) sees the post-cond th and rv!;
      //       this revert could probably be avoided if sstp_tmp_rv/th, pdrv, sstp_save() and etc would be cleverly used
      //       so that th and rv are not directly changed in step_sync (would also fix the problem with chem and src 2 lines above)

      pimpl->should_now_run_async = true;
      pimpl->selected_before_counting = false;
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::step_async(
      const opts_t<real_t> &opts
    ) {
      //sanity checks
      if (!pimpl->should_now_run_async)
        throw std::runtime_error("please call step_sync() before calling step_async() again");

      pimpl->should_now_run_async = false;

      //sanity checks
      if((opts.chem_dsl || opts.chem_dsc || opts.chem_rct) && !pimpl->opts_init.chem_switch) throw std::runtime_error("all chemistry was switched off in opts_init");
      if(opts.coal && !pimpl->opts_init.coal_switch) throw std::runtime_error("all coalescence was switched off in opts_init");
      if(opts.sedi && !pimpl->opts_init.sedi_switch) throw std::runtime_error("all sedimentation was switched off in opts_init");
      if((opts.chem_dsl || opts.chem_dsc || opts.chem_rct) && !pimpl->opts_init.chem_switch) 
        throw std::runtime_error("all chemistry was switched off in opts_init");

      if(opts.coal && !pimpl->opts_init.coal_switch) 
        throw std::runtime_error("all coalescence was switched off in opts_init");

      if(opts.sedi && !pimpl->opts_init.sedi_switch) 
        throw std::runtime_error("all sedimentation was switched off in opts_init");

      if (opts.chem_dsl) 
      { 
        // saving rv to be used as rv_old
        // NOTE: doing it here assumes that gases didn't change since chemistry finished in step_cond
        pimpl->sstp_save_chem();
      }

      // updating Tpr look-up table (includes RH update)
      pimpl->hskpng_Tpr(); 

      // updating terminal velocities
      if (opts.sedi || opts.coal)
        pimpl->hskpng_vterm_all();

      // coalescence
      if (opts.coal) 
      {
        thrust::fill(pimpl->delta_accr20.begin(), pimpl->delta_accr20.end(), real_t(0));
        thrust::fill(pimpl->delta_acnv20.begin(), pimpl->delta_acnv20.end(), real_t(0));
        thrust::fill(pimpl->delta_accr32.begin(), pimpl->delta_accr32.end(), real_t(0));
        thrust::fill(pimpl->delta_acnv32.begin(), pimpl->delta_acnv32.end(), real_t(0));

        for (int step = 0; step < pimpl->opts_init.sstp_coal; ++step) 
        {
          // collide
          pimpl->coal(pimpl->opts_init.dt / pimpl->opts_init.sstp_coal);

          // update invalid vterm 
          if (step + 1 != pimpl->opts_init.sstp_coal)
            pimpl->hskpng_vterm_invalid(); 
        }

        // decrease coalescence timestep
        // done if number of collisions > 1 in const_multi mode
        if(*(pimpl->increase_sstp_coal))
        {
          ++(pimpl->opts_init.sstp_coal);
          *(pimpl->increase_sstp_coal) = false;
        }

        // add accr and acnv to cell statistics, similar to moms_calc + repeated 4 times :)
        {
          const auto n = thrust::reduce_by_key(
            // input - keys
            pimpl->sorted_ijk.begin(), pimpl->sorted_ijk.end(),  
            // input - values
            thrust::make_permutation_iterator(pimpl->delta_accr20.begin(), pimpl->sorted_id.begin()),
            // output - keys
            pimpl->count_ijk.begin(),
            // output - values
            pimpl->count_mom.begin()
          );  
          auto count_n = n.first - pimpl->count_ijk.begin();
          thrust::transform(
            pimpl->count_mom.begin(), pimpl->count_mom.begin() + count_n,
            thrust::make_permutation_iterator(pimpl->accr20.begin(), pimpl->count_ijk.begin()),
            thrust::make_permutation_iterator(pimpl->accr20.begin(), pimpl->count_ijk.begin()),
            thrust::plus<real_t>()
          );
        }
        {
          const auto n = thrust::reduce_by_key(
            // input - keys
            pimpl->sorted_ijk.begin(), pimpl->sorted_ijk.end(),  
            // input - values
            thrust::make_permutation_iterator(pimpl->delta_accr32.begin(), pimpl->sorted_id.begin()),
            // output - keys
            pimpl->count_ijk.begin(),
            // output - values
            pimpl->count_mom.begin()
          );  
          auto count_n = n.first - pimpl->count_ijk.begin();
          thrust::transform(
            pimpl->count_mom.begin(), pimpl->count_mom.begin() + count_n,
            thrust::make_permutation_iterator(pimpl->accr32.begin(), pimpl->count_ijk.begin()),
            thrust::make_permutation_iterator(pimpl->accr32.begin(), pimpl->count_ijk.begin()),
            thrust::plus<real_t>()
          );
        }
        {
          const auto n = thrust::reduce_by_key(
            // input - keys
            pimpl->sorted_ijk.begin(), pimpl->sorted_ijk.end(),  
            // input - values
            thrust::make_permutation_iterator(pimpl->delta_acnv20.begin(), pimpl->sorted_id.begin()),
            // output - keys
            pimpl->count_ijk.begin(),
            // output - values
            pimpl->count_mom.begin()
          );  
          auto count_n = n.first - pimpl->count_ijk.begin();
          thrust::transform(
            pimpl->count_mom.begin(), pimpl->count_mom.begin() + count_n,
            thrust::make_permutation_iterator(pimpl->acnv20.begin(), pimpl->count_ijk.begin()),
            thrust::make_permutation_iterator(pimpl->acnv20.begin(), pimpl->count_ijk.begin()),
            thrust::plus<real_t>()
          );
        }
        {
          const auto n = thrust::reduce_by_key(
            // input - keys
            pimpl->sorted_ijk.begin(), pimpl->sorted_ijk.end(),  
            // input - values
            thrust::make_permutation_iterator(pimpl->delta_acnv32.begin(), pimpl->sorted_id.begin()),
            // output - keys
            pimpl->count_ijk.begin(),
            // output - values
            pimpl->count_mom.begin()
          );  
          auto count_n = n.first - pimpl->count_ijk.begin();
          thrust::transform(
            pimpl->count_mom.begin(), pimpl->count_mom.begin() + count_n,
            thrust::make_permutation_iterator(pimpl->acnv32.begin(), pimpl->count_ijk.begin()),
            thrust::make_permutation_iterator(pimpl->acnv32.begin(), pimpl->count_ijk.begin()),
            thrust::plus<real_t>()
          );
        }
      }

      // advection, it invalidates i,j,k and ijk!
      if (opts.adve) pimpl->adve(); 
      // revert to the desired adve scheme (in case we used eulerian this timestep for halo reasons)
      pimpl->adve_scheme = pimpl->opts_init.adve_scheme;

      // sedimentation has to be done after advection, so that negative z doesnt crash hskpng_ijk in adve
      if (opts.sedi) 
      {
        // advection with terminal velocity
        pimpl->sedi();
      }

      // boundary condition + accumulated rainfall to be returned
      // multi_GPU version invalidates i and k;
      // this has to be done last since i and k will be used by multi_gpu copy to other devices
      // TODO: instead of using i and k define new vectors ?
      // TODO: do this only if we advect/sediment?
      pimpl->bcnd();

      // some stuff to be done at the end of the step.
      // if using more than 1 GPU
      // has to be done after copy 
      if (pimpl->opts_init.dev_count < 2)
        pimpl->step_finalize(opts);

      pimpl->selected_before_counting = false;
    }
  };
};
