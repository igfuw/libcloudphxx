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
        if (!courant_x.is_null()) pimpl->init_e2l(courant_x, &pimpl->courant_x, 1, 0, 0);
        if (!courant_y.is_null()) pimpl->init_e2l(courant_y, &pimpl->courant_y, 0, 1, 0, pimpl->n_x_bfr * pimpl->opts_init.nz);
        if (!courant_z.is_null()) pimpl->init_e2l(courant_z, &pimpl->courant_z, 0, 0, 1, pimpl->n_x_bfr * max(1, pimpl->opts_init.ny));
      }

      // syncing in Eulerian fields (if not null)
      pimpl->sync(th,             pimpl->th);
      pimpl->sync(rv,             pimpl->rv);
      pimpl->sync(courant_x,      pimpl->courant_x);
      pimpl->sync(courant_y,      pimpl->courant_y);
      pimpl->sync(courant_z,      pimpl->courant_z);
      pimpl->sync(rhod,           pimpl->rhod);

      if (pimpl->opts_init.chem_switch){
        for (int i = 0; i < chem_gas_n; ++i){
          pimpl->sync(
            ambient_chem.at((chem_species_t)i), 
            pimpl->ambient_chem[(chem_species_t)i]
          );
        }
      }

      // condensation/evaporation 
      if (opts.cond) 
      {
        for (int step = 0; step < pimpl->opts_init.sstp_cond; ++step) 
        {   
          pimpl->sstp_step(step, !rhod.is_null());
          pimpl->hskpng_Tpr(); 
          pimpl->cond(pimpl->opts_init.dt / pimpl->opts_init.sstp_cond, opts.RH_max); 
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

      // chemistry
      if (opts.chem_dsl or opts.chem_dsc or opts.chem_rct) 
      {
        // calculate new volume of droplets (needed for Henrys law)
        pimpl->chem_vol_ante();
        // set flag for those SD that are big enough to have chemical reactions
        pimpl->chem_flag_ante();

        for (int step = 0; step < pimpl->opts_init.sstp_chem; ++step)
        {
          //dissolving trace gases (Henrys law)
          if (opts.chem_dsl == true)
            pimpl->chem_henry(pimpl->opts_init.dt / pimpl->opts_init.sstp_chem, opts.chem_sys_cls);

          //dissociation
          if (opts.chem_dsc == true)
            pimpl->chem_dissoc();

          //oxidation 
          if (opts.chem_rct == true)
          pimpl->chem_react(pimpl->opts_init.dt / pimpl->opts_init.sstp_chem);
        }

        //save the current drop volume in V_old (to be used in the next step for Henrys law)
        pimpl->chem_vol_post();

        // syncing out // TODO: this is not necesarry in off-line mode (see coupling with DALES)
        for (int i = 0; i < chem_gas_n; ++i)
          pimpl->sync(
            pimpl->ambient_chem[(chem_species_t)i],
            ambient_chem.at((chem_species_t)i)
          );
      }

      pimpl->should_now_run_async = true;
      pimpl->selected_before_counting = false;
    }

    template <typename real_t, backend_t device>
    real_t particles_t<real_t, device>::step_async(
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

      if (opts.cond) 
      { 
        // saving rv to be used as rv_old
        pimpl->sstp_save();
      }

      // updating Tpr look-up table (includes RH update)
      pimpl->hskpng_Tpr(); 

      // advection 
      if (opts.adve) pimpl->adve(); 

      // updating terminal velocities
      if (opts.sedi || opts.coal)
        pimpl->hskpng_vterm_all();

      if (opts.sedi) 
      {
        // advection with terminal velocity
        pimpl->sedi();
      }

      // coalescence
      if (opts.coal) 
      {
        for (int step = 0; step < pimpl->opts_init.sstp_coal; ++step) 
        {
          // collide
          pimpl->coal(pimpl->opts_init.dt / pimpl->opts_init.sstp_coal);

          // update invalid vterm 
          if (step + 1 != pimpl->opts_init.sstp_coal)
            pimpl->hskpng_vterm_invalid(); 
        }
      }

      // boundary condition + accumulated rainfall to be returned
      // multi_GPU version invalidates i and k;
      // this has to be done last since i and k will be used by multi_gpu copy to other devices
      // TODO: instead of using i and k define new vectors ?
      // TODO: do this only if we advect/sediment?
      real_t ret = pimpl->bcnd();

      // some stuff to be done at the end of the step.
      // if using distributed memory
      // has to be done after copy 
      if (!pimpl->dist_mem)
        pimpl->step_finalize(opts);

      pimpl->selected_before_counting = false;

      return ret;
    }
  };
};
