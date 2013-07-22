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
      const ptrdiff_t *strides,
      real_t *rhod_th,
      real_t *rhod_rv,
      real_t *rhod 
    )
    {
std::cerr << "\n\n INIT \n\n";
      // sanity checks
      assert(strides != NULL);
      assert(rhod_th != NULL);
      assert(rhod_rv != NULL);
      assert(rhod != NULL);

      // initialising dry radii
      assert(pimpl->opts.dry_distros.size() == 1); // TODO: handle multiple spectra/kappas
      pimpl->init_dry(pimpl->opts.dry_distros.begin()->second); // TODO: document that n_of_lnrd is expected!

      // initialising particle positions
      pimpl->init_xyz();

      // initialising Euleria-Lagrangian index translator
      pimpl->init_e2l(strides);

      // feeding in Eulerian fields
      pimpl->sync(rhod_th, pimpl->rhod_th);
      pimpl->sync(rhod_rv, pimpl->rhod_rv);
      pimpl->sync(rhod,    pimpl->rhod);

      // calculating derived quantitites
      pimpl->init_Tpr(); // ?

      // 
      pimpl->hskpng_Tpr(); 
      pimpl->hskpng_ijk(); // is it needed here? 
      pimpl->init_wet();
    }
  };
};
