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
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::save_liq_ice_content_before_change() 
    {   
      // --- calc liquid water content before src ---
      hskpng_sort();
      reset_guardp(drw_mom3_gp, tmp_device_real_cell); 
      thrust_device::vector<real_t> &drw_mom3 = drw_mom3_gp->get();

      moms_all(); // TODO: select water only? won't change results, because rw2=0 for non-water particles
      moms_calc(rw2.begin(), real_t(3./2.));

      // drw_mom3 = -rw_mom3 ante change
      if(count_n!=n_cell)  
        thrust::fill(drw_mom3.begin(), drw_mom3.end(), real_t(0.));

      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
        thrust::make_permutation_iterator(drw_mom3.begin(), count_ijk.begin()), // output
        thrust::negate<real_t>()
      );

      if(opts_init.ice_switch)
      {
        reset_guardp(d_ice_mass_gp, tmp_device_real_cell); 
        thrust_device::vector<real_t> &d_ice_mass = d_ice_mass_gp->get();

        moms_gt0(ice_a.begin()); // choose ice particles (ice_a>0)
        moms_calc(thrust::make_transform_iterator(
          thrust::make_zip_iterator(thrust::make_tuple(ice_a.begin(), ice_c.begin(), ice_rho.begin())),
          detail::ice_mass<real_t>()
          ),
          real_t(1)
        );
        nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd ice moment) before deposition");

        // fill with 0s if not all cells have particles
        if(count_n!=n_cell) {
          thrust::fill(d_ice_mass.begin(), d_ice_mass.end(), real_t(0.));
          // thrust::fill(ice_mass.begin(), ice_mass.end(), real_t(0.));
        }

        thrust::transform(
          count_mom.begin(), count_mom.begin() + count_n,
          thrust::make_permutation_iterator(d_ice_mass.begin(), count_ijk.begin()),
          thrust::negate<real_t>()
        );
      }
    }
  };  
};
