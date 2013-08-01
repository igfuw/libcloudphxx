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
    template <typename real_t, int device>
    void particles<real_t, device>::init(
      const arrinfo_t<real_t> rhod_th,
      const arrinfo_t<real_t> rhod_rv,
      const arrinfo_t<real_t> rhod,
      const arrinfo_t<real_t> courant_x, // might be NULL
      const arrinfo_t<real_t> courant_y, // might be NULL
      const arrinfo_t<real_t> courant_z  // might be NULL
    )
    {
std::cerr << "\n\n INIT \n\n";
      // sanity checks
      assert(!rhod_th.is_null());
      assert(!rhod_rv.is_null());
      assert(!rhod.is_null());
      if (pimpl->n_dims > 0) assert(!courant_z.is_null());
      if (pimpl->n_dims > 1) assert(!courant_x.is_null());
      if (pimpl->n_dims > 2) assert(!courant_y.is_null());

      // initialising dry radii
      assert(pimpl->opts.dry_distros.size() == 1); // TODO: handle multiple spectra/kappas
      pimpl->init_dry(pimpl->opts.dry_distros.begin()->second); // TODO: document that n_of_lnrd is expected!

      // initialising particle positions
      pimpl->init_xyz();

      // initialising Eulerian-Lagrandian coupling
      pimpl->init_sync();
std::cerr << "initing rhod_th" << std::endl;
      pimpl->init_e2l(rhod_th, &pimpl->rhod_th);
std::cerr << "initing rhod_rv" << std::endl;
      pimpl->init_e2l(rhod_rv, &pimpl->rhod_rv);
std::cerr << "initing rhod" << std::endl;
      pimpl->init_e2l(rhod,    &pimpl->rhod);
std::cerr << "initing cx" << std::endl;
      if (!courant_x.is_null()) pimpl->init_e2l(courant_x, &pimpl->courant_x, 1, 0, 0);
std::cerr << "initing cy" << std::endl;
      if (!courant_y.is_null()) pimpl->init_e2l(courant_y, &pimpl->courant_y, 0, 1, 0);
std::cerr << "initing cz" << std::endl;
      if (!courant_z.is_null()) pimpl->init_e2l(courant_z, &pimpl->courant_z, 0, 0, 1);

      // feeding in Eulerian fields
      pimpl->sync(rhod_th, pimpl->rhod_th);
      pimpl->sync(rhod_rv, pimpl->rhod_rv);
      pimpl->sync(rhod,    pimpl->rhod);
      if (!courant_x.is_null()) pimpl->sync(courant_x, pimpl->courant_x);
      if (!courant_y.is_null()) pimpl->sync(courant_y, pimpl->courant_y);
      if (!courant_z.is_null()) pimpl->sync(courant_z, pimpl->courant_z);

      // initialising wet radii
      pimpl->init_hskpng(); 
      pimpl->hskpng_Tpr(); 
      pimpl->hskpng_ijk(); 
      pimpl->init_wet();
    }
  };
};
