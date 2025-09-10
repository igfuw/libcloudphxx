// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <thrust/iterator/transform_iterator.h>
#include <libcloudph++/common/detail/toms748.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      // The condition for immersion freezing
      template<class real_t>
      class immersion_freeze_cond
      {
      public:
        BOOST_GPU_ENABLED
        bool operator()(const thrust::tuple<real_t, real_t, real_t> &tpl)  // tpl is a tuple of 3 elements: (T_freeze, ambient T, ambient RH)
        {
          if (thrust::get<0>(tpl) >=  thrust::get<1>(tpl) && thrust::get<2>(tpl) >= real_t(1))  // returns true if T_freeze >= ambient T and ambient RH >= 1
            return true;
          else
            return false;
        };
      };
      // Functor to return rw2 of newly frozen droplets, otherwise 0
      template<class real_t>
      class rw2_if_newfrozen
      {
      public:
        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, int, int> &tpl) const // tpl is a tuple of 3 elements: (rw2, ice, ice_old)
        {
          real_t rw2 = thrust::get<0>(tpl);
          int ice_flag = thrust::get<1>(tpl);
          int ice_flag_old = thrust::get<2>(tpl);
          return (ice_flag==1 && ice_flag_old==0) ? rw2 : real_t(0);
        }
      };

    };

    // Immersion freezing
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::ice_nucl() {

      hskpng_sort();

      // Copy current ice flags
      thrust_device::vector<int> ice_old = ice;

      // Change liquid droplets to ice under the freezing condition
      thrust::replace_if(ice.begin(), ice.end(),                        // Replacing values of ice with 1 if immersion_freeze_cond is satisfied.
        thrust::make_zip_iterator(
          thrust::make_tuple(                                           // Creating a zip iterator to access multiple vectors:
            T_freeze.begin(),                                               // freezing temperature for each droplet
            thrust::make_permutation_iterator(T.begin(), ijk.begin()),      // ambient temperature
            thrust::make_permutation_iterator(RH.begin(), ijk.begin())      // ambient RH
          )
        ),
        detail::immersion_freeze_cond<real_t>(),
        real_t(1)
      );

      // Temporary vector for rw2 of newly frozen droplets
      thrust_device::vector<real_t> rw2_frozen(count_n);

      thrust::transform(
        thrust::make_zip_iterator(thrust::make_tuple( // first input
          rw2.begin(),    // droplet radius squared
          ice.begin(),    // ice flag
          ice_old.begin() // old ice flag
        )),
        thrust::make_zip_iterator(thrust::make_tuple( // last input
          rw2.end(),
          ice.end(),
          ice_old.end()
        )),
        rw2_frozen.begin(),                 // output
        detail::rw2_if_newfrozen<real_t>()   // functor for r2 of newly frozen droplets
      );

      // Compute per-cell 3rd moment of newly frozen droplets (sum of n*r^3). It is stored in count_mom
      moms_all();
      moms_calc(rw2_frozen.begin(), count_n, real_t(1.5), true);
      nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd wet moment) of newly frozen droplets");

      // Copy to frozen_mom3 array
      thrust_device::vector<real_t> &frozen_mom3(tmp_device_real_cell);
      if (count_n != n_cell)
        thrust::fill(frozen_mom3.begin(), frozen_mom3.end(), real_t(0.));

      thrust::transform(
          count_mom.begin(), count_mom.begin() + count_n,
          thrust::make_permutation_iterator(frozen_mom3.begin(), count_ijk.begin()),
          thrust::identity<real_t>()  // just copy values
      );

      // update th according to the frozen volume per cell
      update_th_freezing(frozen_mom3);
    }

  }
}
