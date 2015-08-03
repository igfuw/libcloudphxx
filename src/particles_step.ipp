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
      const arrinfo_t<real_t> courant_x, // defaults to NULL-NULL pair (e.g. kinematic model)
      const arrinfo_t<real_t> courant_y, // defaults to NULL-NULL pair (e.g. kinematic model)
      const arrinfo_t<real_t> courant_z, // defaults to NULL-NULL pair (e.g. kinematic model)
      const arrinfo_t<real_t> rhod       // defaults to NULL-NULL pair (e.g. kinematic or boussinesq model)
    )
    {
      if (!pimpl->init_called)
        throw std::runtime_error("please call init() before calling step_sync()");
      if (pimpl->should_now_run_async)
        throw std::runtime_error("please call step_async() before calling step_sync() again");

      if (pimpl->l2e[&pimpl->courant_x].size() == 0) // TODO: y, z,...
      {
        // TODO: copy-pasted from init
	if (!courant_x.is_null()) pimpl->init_e2l(courant_x, &pimpl->courant_x, 1, 0, 0); 
	if (!courant_y.is_null()) pimpl->init_e2l(courant_y, &pimpl->courant_y, 0, 1, 0); 
	if (!courant_z.is_null()) pimpl->init_e2l(courant_z, &pimpl->courant_z, 0, 0, 1);
      }

      // syncing in Eulerian fields (if not null)
      pimpl->sync(th,             pimpl->th);
      pimpl->sync(rv,             pimpl->rv);
      pimpl->sync(courant_x,      pimpl->courant_x);
      pimpl->sync(courant_y,      pimpl->courant_y);
      pimpl->sync(courant_z,      pimpl->courant_z);
      pimpl->sync(rhod,           pimpl->rhod);

      // recycling out-of-domain/invalidated particles 
      // (doing it here and not in async reduces the need for a second sort before diagnostics,
      // but also unneccesarily holds dyncore execution for a bit longer)
      thrust_size_t n_rcyc = 0;//pimpl->rcyc();

      // updating particle->cell look-up table
      // (before advection and sedimentation so that their order does not matter,
      if (opts.adve || opts.sedi || n_rcyc)
      {
        pimpl->hskpng_ijk();
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

        // syncing out // TODO: this is not necesarry in off-line mode (see coupling with DALES)
        pimpl->sync(pimpl->th, th);
        pimpl->sync(pimpl->rv, rv);
      }

      pimpl->should_now_run_async = true;
      pimpl->selected_before_counting = false;
    }

    template <typename real_t, backend_t device>
    real_t particles_t<real_t, device>::step_async(
      const opts_t<real_t> &opts
    ) {
      if (!pimpl->should_now_run_async)
        throw std::runtime_error("please call step_sync() before calling step_async() again");

      //sanity checks
      if((opts.chem_dsl || opts.chem_dsc || opts.chem_rct) && !pimpl->opts_init.chem_switch) throw std::runtime_error("all chemistry was switched off in opts_init");
      if(opts.coal && !pimpl->opts_init.coal_switch) throw std::runtime_error("all coalescence was switched off in opts_init");
      if(opts.sedi && !pimpl->opts_init.sedi_switch) throw std::runtime_error("all sedimentation was switched off in opts_init");

      if (opts.cond) 
      { 
        // saving rv to be used as rv_old
        pimpl->sstp_save();
      }

      // updating Tpr look-up table (includes RH update)
      pimpl->hskpng_Tpr(); 

      // advection 
      if (opts.adve) pimpl->adve(); 

      // boundary condition + accumulated rainfall to be returned
      real_t ret = pimpl->bcnd();

      // updating terminal velocities
      if (opts.sedi || opts.coal)
        pimpl->hskpng_vterm_all();

      if (opts.sedi) 
      {
        // advection with terminal velocity
        pimpl->sedi();
      }

      // chemistry
      if (opts.chem_dsl or opts.chem_dsc or opts.chem_rct) 
      {
        for (int step = 0; step < pimpl->opts_init.sstp_chem; ++step) 
          pimpl->chem(pimpl->opts_init.dt / pimpl->opts_init.sstp_chem, opts.chem_gas, 
                      opts.chem_dsl, opts.chem_dsc, opts.chem_rct
                     );
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

      pimpl->should_now_run_async = false;
      pimpl->selected_before_counting = false;

      // aerosol source
      if (opts.src) 
      {
        if (pimpl->opts_init.src_switch == false) throw std::runtime_error("aerosol source was switched off in opts_init");

        pimpl->src();
      }

      return ret;
    }
  };
};
