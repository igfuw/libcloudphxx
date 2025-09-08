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
      // Functor to compute r^3 only for newly frozen particles
      template<class real_t>
      class compute_r3_if_frozen
      {
      public:
        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, int, int> &tpl) const // tpl is a tuple of 3 elements: (r2, ice, ice_old)
        {
          real_t r2 = thrust::get<0>(tpl);
          int ice_flag = thrust::get<1>(tpl);
          int ice_flag_old = thrust::get<2>(tpl);
          return (ice_flag==1 && ice_flag_old==0) ? pow(r2, real_t(1.5)) : real_t(0);
        }
      };

    };

    // Immersion freezing
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::ice_nucl() {

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

      // Temporary vector for 3rd moment contribution of each droplet
      thrust_device::vector<real_t> mom3(count_n);

      // Compute r^3 only for newly frozen droplets
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
        mom3.begin(),                 // output
        detail::compute_r3_if_frozen<real_t>()   // functor for computing r3 of newly frozen droplets
      );

      // Reuse a temporary device vector for cell-wise 3rd moment
      thrust_device::vector<real_t> &dri(tmp_device_real_cell);
      thrust::fill(dri.begin(), dri.end(), real_t(0));  // reset to 0

      // add contributions to cell-wise 3rd moment
      thrust::transform(
        mom3.begin(), mom3.end(),                       // input
        thrust::make_permutation_iterator(dri.begin(), count_ijk.begin()), // output per cell
        thrust::make_permutation_iterator(dri.begin(), count_ijk.begin()),
        thrust::plus<real_t>()
      );

      // update th according to changes in ri
      update_th_freezing(dri);
    }

  }
}
