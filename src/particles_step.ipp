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
      const arrinfo_t<real_t> rhod_courant_x, // defaults to NULL-NULL pair (e.g. kinematic model)
      const arrinfo_t<real_t> rhod_courant_y, // defaults to NULL-NULL pair (e.g. kinematic model)
      const arrinfo_t<real_t> rhod_courant_z, // defaults to NULL-NULL pair (e.g. kinematic model)
      const arrinfo_t<real_t> rhod       // defaults to NULL-NULL pair (e.g. kinematic or boussinesq model)
    )
    {
      assert(!pimpl->should_now_run_async);

      // syncing in Eulerian fields (if not null)
      pimpl->sync(th,             pimpl->th);
      pimpl->sync(rv,             pimpl->rv);
      pimpl->sync(rhod_courant_x, pimpl->rhod_courant_x);
      pimpl->sync(rhod_courant_y, pimpl->rhod_courant_y);
      pimpl->sync(rhod_courant_z, pimpl->rhod_courant_z);
      pimpl->sync(rhod,           pimpl->rhod);

      // recycling out-of-domain/invalidated particles 
      // (doing it here and not in async reduces the need for a second sort before diagnostics,
      // but also unneccesarily holds dyncore execution for a bit longer)
      pimpl->rcyc(); 

      // updating particle->cell look-up table
      // (before advection and sedimentation so that their order does not matter,
      if (opts.adve || opts.sedi)
      {
        pimpl->hskpng_ijk(); // TODO: but rcyc() above could also have changed ijk!
      }

      // updating Tpr look-up table (includes RH update)
      pimpl->hskpng_Tpr(); 

      // condensation/evaporation 
      if (opts.cond) 
      {
        for (int step = 0; step < opts.sstp_cond; ++step) 
        {   
          pimpl->cond(pimpl->opts_init.dt / opts.sstp_cond, opts.RH_max); 
          pimpl->hskpng_Tpr(); 
        } 

        // syncing out 
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
      assert(pimpl->should_now_run_async);

      // accumulated rainfall to be returned
      real_t ret;

      // advection 
      if (opts.adve) pimpl->adve(); 

      // boundary condition
      ret = pimpl->bcnd();

      // updating terminal velocities
      if (opts.sedi || opts.coal) 
        pimpl->hskpng_vterm_all();

      if (opts.sedi) 
      {
        // advection with terminal velocity
        pimpl->sedi();
      }

      // chemistry
      if (opts.chem) 
      {
        for (int step = 0; step < opts.sstp_chem; ++step) 
          pimpl->chem(pimpl->opts_init.dt / opts.sstp_chem, opts.chem_gas);
      }

      // coalescence (before diagnostics -> one sort less)
      if (opts.coal) 
      {
        for (int step = 0; step < opts.sstp_coal; ++step) 
        {
          // collide
          pimpl->coal(pimpl->opts_init.dt / opts.sstp_coal);

          // update invalid vterm 
          if (step + 1 != opts.sstp_coal)
            pimpl->hskpng_vterm_invalid(); 
        }
      }

      pimpl->should_now_run_async = false;
      pimpl->selected_before_counting = false;

      return ret;
    }
  };
};
