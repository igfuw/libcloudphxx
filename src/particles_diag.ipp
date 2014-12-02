// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

namespace libcloudphxx
{
  namespace lgrngn
  {
    // records super-droplet concentration per grid cell
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_sd_conc()
    {
      // common code with coalescence, hence separated into a method
      pimpl->hskpng_count(); 

      // n_t -> real_t cast
      thrust::copy(
        pimpl->count_num.begin(), 
        pimpl->count_num.end(), 
        pimpl->count_mom.begin()
      );
    }

    // selected all particles
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_all()
    {
      pimpl->moms_all();
    }

    // selects particles with (r_d >= r_min && r_d < r_max)
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_dry_rng(const real_t &r_min, const real_t &r_max)
    {
      pimpl->moms_rng(pow(r_min, 3), pow(r_max, 3), pimpl->rd3.begin());
    }

    // selects particles with (r_w >= r_min && r_w < r_max)
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_wet_rng(const real_t &r_min, const real_t &r_max)
    {
      pimpl->moms_rng(pow(r_min, 2), pow(r_max, 2), pimpl->rw2.begin());
    }

    // computes n-th moment of the dry spectrum for the selected particles
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_dry_mom(const int &n)
    {
      pimpl->moms_calc(pimpl->rd3.begin(), n/3.);
    }

    // computes n-th moment of the wet spectrum for the selected particles
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_wet_mom(const int &n)
    {
      pimpl->moms_calc(pimpl->rw2.begin(), n/2.);
    }

    // computes mean chemical properties for the selected particles
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::diag_chem(const enum chem_species_t &c)
    {
      pimpl->moms_calc(pimpl->chem_bgn[c], 1);
    }
  };
};
