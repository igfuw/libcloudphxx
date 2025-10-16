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
      const arrinfo_t<real_t> diss_rate, // defaults to NULL-NULL pair
      std::map<enum chem_species_t, arrinfo_t<real_t> > ambient_chem
    )
    {
      sync_in(th, rv, rhod, courant_x, courant_y, courant_z, diss_rate, ambient_chem);
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
      const arrinfo_t<real_t> diss_rate, // defaults to NULL-NULL pair
      std::map<enum chem_species_t, arrinfo_t<real_t> > ambient_chem
    )
    {
      // sanity checks
      if (!pimpl->init_called)
        throw std::runtime_error("libcloudph++: please call init() before calling step_sync()");
      if (pimpl->should_now_run_async)
        throw std::runtime_error("libcloudph++: please call step_async() before calling step_sync() again");

      if (th.is_null() || rv.is_null())
        throw std::runtime_error("libcloudph++: passing th and rv is mandatory");

 // <TODO> - code duplicated from init() !
      if (!courant_x.is_null() || !courant_y.is_null() || !courant_z.is_null())
      {
	if (pimpl->n_dims == 0)
	  throw std::runtime_error("libcloudph++: Courant numbers passed in 0D setup");

	if (pimpl->n_dims == 1 && (courant_x.is_null() || !courant_y.is_null() || !courant_z.is_null()))
	  throw std::runtime_error("libcloudph++: Only X Courant number allowed in 1D setup");

	if (pimpl->n_dims == 2 && (courant_x.is_null() || !courant_y.is_null() || courant_z.is_null()))
	  throw std::runtime_error("libcloudph++: Only X and Z Courant numbers allowed in 2D setup");

	if (pimpl->n_dims == 3 && (courant_x.is_null() || courant_y.is_null() || courant_z.is_null()))
	  throw std::runtime_error("libcloudph++: All XYZ Courant number components required in 3D setup");
      }

      if (pimpl->opts_init.chem_switch && ambient_chem.size() != chem_gas_n)
        throw std::runtime_error("libcloudph++: chemistry was not switched off and ambient_chem is empty");

      if (!pimpl->opts_init.chem_switch && ambient_chem.size() != 0)
        throw std::runtime_error("libcloudph++: chemistry was switched off and ambient_chem is not empty");

      if ( (pimpl->opts_init.turb_adve_switch || pimpl->opts_init.turb_cond_switch || pimpl->opts_init.turb_coal_switch)  && diss_rate.is_null())
        throw std::runtime_error("libcloudph++: turbulent advection, coalescence and condesation are not switched off and diss_rate is empty");

      if ( !(pimpl->opts_init.turb_adve_switch || pimpl->opts_init.turb_cond_switch || pimpl->opts_init.turb_coal_switch)  && !diss_rate.is_null())
        throw std::runtime_error("libcloudph++: turbulent advection, coalescence and condesation are switched off and diss_rate is not empty");
// </TODO>

      if (pimpl->l2e[&pimpl->courant_x].size() == 0) // TODO: y, z,...
      {
        // TODO: copy-pasted from init
        if (!courant_x.is_null()) pimpl->init_e2l(courant_x, &pimpl->courant_x, 1, 0, 0, - pimpl->halo_x);
        if (!courant_y.is_null()) pimpl->init_e2l(courant_y, &pimpl->courant_y, 0, 1, 0, pimpl->n_x_bfr * pimpl->opts_init.nz - pimpl->halo_y);
        if (!courant_z.is_null()) pimpl->init_e2l(courant_z, &pimpl->courant_z, 0, 0, 1, pimpl->n_x_bfr * std::max(1, pimpl->opts_init.ny) - pimpl->halo_z);
        /*
        if (!courant_x.is_null()) pimpl->init_e2l(courant_x, &pimpl->courant_x, 1, 0, 0);
        if (!courant_y.is_null()) pimpl->init_e2l(courant_y, &pimpl->courant_y, 0, 1, 0, pimpl->n_x_bfr * pimpl->opts_init.nz );
        if (!courant_z.is_null()) pimpl->init_e2l(courant_z, &pimpl->courant_z, 0, 0, 1, pimpl->n_x_bfr * std::max(1, pimpl->opts_init.ny) );
        */
      }

      if (pimpl->l2e[&pimpl->diss_rate].size() == 0)
        if (!diss_rate.is_null()) pimpl->init_e2l(diss_rate, &pimpl->diss_rate);

      // did rhod change
      pimpl->var_rho = !rhod.is_null();

      // syncing in Eulerian fields (if not null)
      pimpl->sync(th,             pimpl->th);
      pimpl->sync(rv,             pimpl->rv);
      pimpl->sync(diss_rate,      pimpl->diss_rate);
      pimpl->sync(rhod,           pimpl->rhod);
      pimpl->sync(courant_x,      pimpl->courant_x);
      pimpl->sync(courant_y,      pimpl->courant_y);
      pimpl->sync(courant_z,      pimpl->courant_z);

      // fill in mpi courant halos
      pimpl->xchng_courants();

      nancheck(pimpl->th, " th after sync-in");
      nancheck(pimpl->rv, " rv after sync-in");
      nancheck(pimpl->courant_x, " courant_x after sync-in");
      nancheck(pimpl->courant_y, " courant_y after sync-in");
      nancheck(pimpl->courant_z, " courant_z after sync-in");
      nancheck(pimpl->diss_rate, " diss_rate after sync-in");
      nancheck(pimpl->rhod, " rhod after sync-in");
      if(pimpl->opts_init.turb_adve_switch || pimpl->opts_init.turb_cond_switch || pimpl->opts_init.turb_coal_switch)
        {nancheck(pimpl->diss_rate, " diss_rate after sync-in");}

      assert(*thrust::min_element(pimpl->rv.begin(), pimpl->rv.end()) >= 0);
      assert(*thrust::min_element(pimpl->th.begin(), pimpl->th.end()) >= 0);
      assert(*thrust::min_element(pimpl->rhod.begin(), pimpl->rhod.end()) >= 0);
      if(pimpl->opts_init.turb_adve_switch || pimpl->opts_init.turb_cond_switch || pimpl->opts_init.turb_coal_switch)
        {assert(*thrust::min_element(pimpl->diss_rate.begin(), pimpl->diss_rate.end()) >= 0);}

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
          nancheck_range(pimpl->ambient_chem[(chem_species_t)i].begin(), pimpl->ambient_chem[(chem_species_t)i].end(), "ambient_chem after sync in");
          assert(*thrust::min_element(pimpl->ambient_chem[(chem_species_t)i].begin(), pimpl->ambient_chem[(chem_species_t)i].end()) >= 0);
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
        throw std::runtime_error("libcloudph++: please call sync_in() before calling step_cond()");

      if(opts.turb_cond && !pimpl->opts_init.turb_cond_switch) 
        throw std::runtime_error("libcloudph++: turb_cond_swtich=False, but turb_cond==True");

      pimpl->should_now_run_cond = false;

      // dt defined in opts_init can be overriden by dt in opts
      pimpl->adjust_timesteps(opts.dt);

      if (pimpl->opts_init.diag_incloud_time)
        pimpl->update_incloud_time(pimpl->dt);

      // condensation/evaporation 
      if (opts.cond) 
      {
        // calculate mean free path needed for molecular correction
        // NOTE: this should be don per each substep, but right now there is no logic
        //       that would make it easy to do in exact (per-cell) substepping
        pimpl->hskpng_mfp(); 

        if(pimpl->opts_init.exact_sstp_cond && pimpl->sstp_cond > 1)
        // apply substeps per-particle logic
        {
          for (int step = 0; step < pimpl->sstp_cond; ++step) 
          {   
            pimpl->sstp_step_exact(step);
            if(opts.turb_cond)
              pimpl->sstp_step_ssp(pimpl->dt / pimpl->sstp_cond);
            if (opts.ice_nucl)
              pimpl->ice_nucl_melt(pimpl->dt / pimpl->sstp_cond, opts.time_dep_ice_nucl);
            pimpl->cond_sstp(pimpl->dt / pimpl->sstp_cond, opts.RH_max, opts.turb_cond);
          } 
          // copy sstp_tmp_rv and th to rv and th
          pimpl->update_state(pimpl->rv, pimpl->sstp_tmp_rv);
          pimpl->update_state(pimpl->th, pimpl->sstp_tmp_th);
        }
        else
        // apply per-cell sstp logic
        {
          for (int step = 0; step < pimpl->sstp_cond; ++step) 
          {   
            pimpl->sstp_step(step);
            if(opts.turb_cond)
              pimpl->sstp_step_ssp(pimpl->dt / pimpl->sstp_cond);
            pimpl->hskpng_Tpr();
            if (opts.ice_nucl)
              pimpl->ice_nucl_melt(pimpl->dt / pimpl->sstp_cond, opts.time_dep_ice_nucl);
            pimpl->cond(pimpl->dt / pimpl->sstp_cond, opts.RH_max, opts.turb_cond);
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
        for (int step = 0; step < pimpl->sstp_chem; ++step) 
        {   
          // calculate new volume of droplets (needed for chemistry)
          pimpl->chem_vol_ante();
          // set flag for those SD that are big enough to have chemical reactions
          pimpl->chem_flag_ante();
          // NOTE: volume and flag are stored in temporary arrays (tmp_device_real_part and tmp_device_n_part), so these arrays should not be changed until chemistry is done

          if (opts.chem_dsl)
          {
            //adjust trace gases to substepping
            pimpl->sstp_step_chem(step);

            //dissolving trace gases (Henrys law)
            pimpl->chem_henry(pimpl->dt / pimpl->sstp_chem);

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
            pimpl->chem_react(pimpl->dt / pimpl->sstp_chem);

            //cleanup - TODO think of something better
            pimpl->chem_cleanup();
          }
        }
      }

      if(opts.cond) // || (opts.src && pimpl->src_stp_ctr % pimpl->opts_init.supstp_src == 0) || (opts.rlx && pimpl->rlx_stp_ctr % pimpl->opts_init.supstp_rlx == 0))
      {
        // syncing out // TODO: this is not necesarry in off-line mode (see coupling with DALES)
        pimpl->sync(pimpl->th, th);
        pimpl->sync(pimpl->rv, rv);
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
        throw std::runtime_error("libcloudph++: please call step_sync() before calling step_async() again");

      pimpl->should_now_run_async = false;

      //sanity checks
      if((opts.chem_dsl || opts.chem_dsc || opts.chem_rct) && !pimpl->opts_init.chem_switch) throw std::runtime_error("libcloudph++: all chemistry was switched off in opts_init");
      if(opts.coal && !pimpl->opts_init.coal_switch) throw std::runtime_error("libcloudph++: coalescence was switched off in opts_init");
      if(opts.sedi && !pimpl->opts_init.sedi_switch) throw std::runtime_error("libcloudph++: sedimentation was switched off in opts_init");
      if(opts.subs && !pimpl->opts_init.subs_switch) throw std::runtime_error("libcloudph++: subsidence was switched off in opts_init");
      if((opts.chem_dsl || opts.chem_dsc || opts.chem_rct) && !pimpl->opts_init.chem_switch) 
        throw std::runtime_error("libcloudph++: all chemistry was switched off in opts_init");

      if(opts.turb_adve && !pimpl->opts_init.turb_adve_switch) 
        throw std::runtime_error("libcloudph++: turb_adve_switch=False, but turb_adve==True");

      if(opts.turb_adve && pimpl->n_dims==0) 
        throw std::runtime_error("libcloudph++: turbulent advection does not work in 0D");

      // dt defined in opts_init can be overriden by dt in opts
      pimpl->adjust_timesteps(opts.dt);

      if (opts.chem_dsl) 
      { 
        // saving rv to be used as rv_old
        // NOTE: doing it here assumes that gases didn't change since chemistry finished in step_cond
        pimpl->sstp_save_chem();
      }

      // updating Tpr look-up table (includes RH update)
      pimpl->hskpng_Tpr(); 

      // updating terminal velocities
      if (opts.sedi || opts.coal || opts.cond)
        pimpl->hskpng_vterm_all();

      // coalescence
      if (opts.coal) 
      {
        for (int step = 0; step < pimpl->sstp_coal; ++step) 
        {
          // collide
          pimpl->coal(pimpl->dt / pimpl->sstp_coal, opts.turb_coal);

          // update invalid vterm 
          if (step + 1 != pimpl->sstp_coal)
            pimpl->hskpng_vterm_invalid(); 
        }

        // decrease coalescence timestep
        // done if number of collisions > 1 in const_multi mode
        if(*(pimpl->increase_sstp_coal))
        {
          ++(pimpl->opts_init.sstp_coal);
          *(pimpl->increase_sstp_coal) = false;
        }
      }

      if (opts.turb_adve || opts.turb_cond)
      {
        // calc tke (diss_rate now holds TKE, not dissipation rate! Hence this must be done after coal, which requires diss rate)
        pimpl->hskpng_tke();
      }
      if (opts.turb_adve)
      {
        // calc turbulent perturbation of velocity
        pimpl->hskpng_turb_vel(pimpl->dt);
      }
      else if (opts.turb_cond)
      {
        // calc turbulent perturbation only of vertical velocity
        pimpl->hskpng_turb_vel(pimpl->dt, true);
      }

      if(opts.turb_cond)
      {
        // calculate the time derivatie of the turbulent supersaturation perturbation; applied in the next step during condensation substepping - is the delay a problem?
        pimpl->hskpng_turb_dot_ss(); 
      }

      // advection, it invalidates i,j,k and ijk!
      if (opts.adve) pimpl->adve(); 
      // revert to the desired adve scheme (in case we used eulerian this timestep for halo reasons)
      pimpl->adve_scheme = pimpl->opts_init.adve_scheme;

      // apply turbulent perturbation of velocity, TODO: add it to advection velocity (turb_vel_calc would need to be called couple times in the pred-corr advection + diss_rate would need a halo)
      if (opts.turb_adve) pimpl->turb_adve(pimpl->dt);

      // sedimentation/subsidence has to be done after advection, so that negative z doesnt crash hskpng_ijk in adve
      if (opts.sedi) 
      {
        // advection with terminal velocity, TODO: add it to the advection velocity (makes a difference for predictor-corrector)
        pimpl->sedi(pimpl->dt);
      }
      if (opts.subs) 
      {
        // advection with subsidence velocity, TODO: add it to the advection velocity (makes a difference for predictor-corrector)
        pimpl->subs(pimpl->dt);
      }

      // NOTE: source and relax should affect th and rv (because we add humidifed aerosols), but these changes are minimal and we neglect them.
      //       otherwise src and rlx should be don in step_sync

      // aerosol source, in sync since it changes th/rv
      if (opts.src && !(pimpl->opts_init.src_x0 == 0 && pimpl->opts_init.src_x1 == 0)) // src_x0=0 and src_x1=0 is a way of disabling source in some domains in distmem simulations
      {
        // sanity check
        if (pimpl->opts_init.src_type == src_t::off) throw std::runtime_error("libcloudph++: aerosol source was switched off in opts_init");

        // introduce new particles
        pimpl->src(opts.src_dry_distros, opts.src_dry_sizes);
      }

      // aerosol relaxation, in sync since it changes th/rv
      // TODO: more sanity checks for rlx! 3D only, values of rlx_bins etc. check that appa ranges are exclusive, zmin<zmax, kpamin<kpamax,...
      if (opts.rlx)
      {
        // sanity check
        if (pimpl->opts_init.rlx_switch == false) throw std::runtime_error("libcloudph++: aerosol relaxation was switched off in opts_init");

        // introduce new particles with the given time interval
        if(pimpl->rlx_stp_ctr % pimpl->opts_init.supstp_rlx == 0) 
        {
          pimpl->rlx(pimpl->opts_init.supstp_rlx * pimpl->dt);
        }
      }

      // update the step counter since src/rlx was turned on
      if (opts.src) ++pimpl->src_stp_ctr;
      else pimpl->src_stp_ctr = 0; //reset the counter if source was turned off
      if (opts.rlx) ++pimpl->rlx_stp_ctr;
      else pimpl->rlx_stp_ctr = 0; //reset the counter if source was turned off

      // boundary condition + accumulated rainfall to be returned
      // distmem version overwrites i and tmp_device_size_part
      // and they both need to be unchanged untill distmem copies
      pimpl->bcnd();
      
      // copy advected SDs using asynchronous MPI;
      if (opts.adve || opts.turb_adve)
        pimpl->mpi_exchange();

      // stuff has to be done after distmem copy 
      // if it is a spawn of multi_CUDA, multi_CUDA will handle finalize
      if(!pimpl->opts_init.dev_count)
        pimpl->post_copy(opts);

      pimpl->selected_before_counting = false;
    }
  };
};
