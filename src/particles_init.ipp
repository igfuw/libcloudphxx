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
      const arrinfo_t<real_t> rhod_courant_x, // might be NULL
      const arrinfo_t<real_t> rhod_courant_y, // might be NULL
      const arrinfo_t<real_t> rhod_courant_z  // might be NULL
    )
    {
      // sanity checks
      assert(!th.is_null());
      assert(!rv.is_null());
      assert(!rhod.is_null());
      if (pimpl->n_dims > 0) assert(!rhod_courant_z.is_null());
      if (pimpl->n_dims > 1) assert(!rhod_courant_x.is_null());
      if (pimpl->n_dims > 2) assert(!rhod_courant_y.is_null());

      // initialising Eulerian-Lagrandian coupling
      pimpl->init_sync();
      pimpl->init_e2l(th,   &pimpl->th);
      pimpl->init_e2l(rv,   &pimpl->rv);
      pimpl->init_e2l(rhod, &pimpl->rhod);

      if (!rhod_courant_x.is_null()) pimpl->init_e2l(rhod_courant_x, &pimpl->rhod_courant_x, 1, 0, 0);
      if (!rhod_courant_y.is_null()) pimpl->init_e2l(rhod_courant_y, &pimpl->rhod_courant_y, 0, 1, 0);
      if (!rhod_courant_z.is_null()) pimpl->init_e2l(rhod_courant_z, &pimpl->rhod_courant_z, 0, 0, 1);

      // feeding in Eulerian fields
      pimpl->sync(th,   pimpl->th);
      pimpl->sync(rv,   pimpl->rv);
      pimpl->sync(rhod, pimpl->rhod);
      if (!rhod_courant_x.is_null()) pimpl->sync(rhod_courant_x, pimpl->rhod_courant_x);
      if (!rhod_courant_y.is_null()) pimpl->sync(rhod_courant_y, pimpl->rhod_courant_y);
      if (!rhod_courant_z.is_null()) pimpl->sync(rhod_courant_z, pimpl->rhod_courant_z);

      // initialising particle positions
      pimpl->init_xyz();

      // initialising helper data for advection (Arakawa-C grid neighbours' indices)
      pimpl->init_grid();

      // initialising housekeeping data (incl. ijk)
      pimpl->init_hskpng(); 
      pimpl->hskpng_Tpr(); 
      pimpl->hskpng_ijk(); 

      // initialising dry radii (needs positions, ijk and rhod)
      assert(pimpl->opts_init.dry_distros.size() == 1); // TODO: handle multiple spectra/kappas
      pimpl->init_dry(
        pimpl->opts_init.dry_distros.begin()->first,
        pimpl->opts_init.dry_distros.begin()->second 
      ); // TODO: document that n_of_lnrd_stp is expected!

      // initialising wet radii
      pimpl->init_wet();
 
      // initialising chem stuff
      pimpl->init_chem(); // TODO: only if chem enabled?
    }
  };
};
