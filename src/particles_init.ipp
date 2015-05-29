// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief initialisation routine for super droplets
  */

namespace libcloudphxx
{
  namespace lgrngn
  {
    // init
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::init(
      const arrinfo_t<real_t> th,
      const arrinfo_t<real_t> rv,
      const arrinfo_t<real_t> rhod,
      const arrinfo_t<real_t> courant_x, // might be NULL
      const arrinfo_t<real_t> courant_y, // might be NULL
      const arrinfo_t<real_t> courant_z  // might be NULL
    )
    {
      // sanity checks
      assert(!th.is_null());
      assert(!rv.is_null());
      assert(!rhod.is_null());

      // initialising Eulerian-Lagrandian coupling
      pimpl->init_sync();
      pimpl->init_e2l(th,   &pimpl->th);
      pimpl->init_e2l(rv,   &pimpl->rv);
      pimpl->init_e2l(rhod, &pimpl->rhod);

      if (!courant_x.is_null()) pimpl->init_e2l(courant_x, &pimpl->courant_x, 1, 0, 0);
      if (!courant_y.is_null()) pimpl->init_e2l(courant_y, &pimpl->courant_y, 0, 1, 0);
      if (!courant_z.is_null()) pimpl->init_e2l(courant_z, &pimpl->courant_z, 0, 0, 1);

      // feeding in Eulerian fields
      pimpl->sync(th,   pimpl->th);
      pimpl->sync(rv,   pimpl->rv);
      pimpl->sync(rhod, pimpl->rhod);
      if (!courant_x.is_null()) pimpl->sync(courant_x, pimpl->courant_x);
      if (!courant_y.is_null()) pimpl->sync(courant_y, pimpl->courant_y);
      if (!courant_z.is_null()) pimpl->sync(courant_z, pimpl->courant_z);

      // initialising housekeeping data of the size of number of cells
      pimpl->init_hskpng_ncell(); 

      // initialising dry radii (needs rhod) with constant multiplicity method (DSMC-like), done before initialization of positions
      if(pimpl->opts_init.sd_const_multi > 0)
      {
        assert(pimpl->opts_init.dry_distros.size() == 1); // TODO: handle multiple spectra/kappas
        pimpl->init_dry_const_multi(
          pimpl->opts_init.dry_distros.begin()->first,
          pimpl->opts_init.dry_distros.begin()->second 
        ); // TODO: document that n_of_lnrd_stp is expected!
      }

      // initialising housekeeping data of the size of number of parts (could have been changed by init_dry_const_multi)
      pimpl->init_hskpng_npart(); 

printf("xyz\n");
      // initialising particle positions
      pimpl->init_xyz();

printf("grid\n");
      // initialising helper data for advection (Arakawa-C grid neighbours' indices)
      pimpl->init_grid();

printf("Tpr\n");
      // initialising additional housekeeping data (incl. ijk)
      pimpl->hskpng_Tpr(); 
printf("ijk\n");
      pimpl->hskpng_ijk(); 

      // initialising dry radii (needs positions, ijk and rhod)
      if(pimpl->opts_init.sd_conc > 0)
      {
        assert(pimpl->opts_init.dry_distros.size() == 1); // TODO: handle multiple spectra/kappas
        pimpl->init_dry(
          pimpl->opts_init.dry_distros.begin()->first,
          pimpl->opts_init.dry_distros.begin()->second 
        ); // TODO: document that n_of_lnrd_stp is expected!
      }

printf("wet\n");
      // initialising wet radii
      pimpl->init_wet();

      // initialising chem stuff
      if(pimpl->opts_init.chem_switch) pimpl->init_chem();

      pimpl->init_sstp();

      //initialising collision kernel
      pimpl->init_kernel();
    }
  };
};
