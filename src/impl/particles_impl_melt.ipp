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
      // The condition for melting
      template<class real_t>
      class melting_cond
      {
      public:
        BOOST_GPU_ENABLED
        bool operator()(const real_t &T) const
        {
          return (T > real_t(273.15));
        };
      };

      // Functor to return rw2 of newly melted particles, otherwise 0
      template<class real_t>
      class rw2_if_newmelted
      {
      public:
        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, int, int> &tpl) const // tpl is a tuple of 3 elements: (rw2, ice, ice_old)
        {
          real_t rw2 = thrust::get<0>(tpl);
          int ice_flag = thrust::get<1>(tpl);
          int ice_flag_old = thrust::get<2>(tpl);
          return (ice_flag==0 && ice_flag_old==1) ? rw2 : real_t(0);
        }
      };

    };

    // Melting
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::melt() {

      hskpng_sort();

      // Copy current ice flags
      thrust_device::vector<int> ice_old = ice;

      // Change ice to liquid under the melting condition
      thrust::replace_if(ice.begin(), ice.end(),                    // Replacing values of ice with 0 if melting_cond is satisfied.
        thrust::make_permutation_iterator(T.begin(), ijk.begin()),  // ambient temperature for each particle
        detail::melting_cond<real_t>(),
        real_t(0)
      );

      // Temporary vector for rw2 of newly melted particles
      thrust_device::vector<real_t> rw2_melted(count_n);

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
        rw2_melted.begin(),                 // output
        detail::rw2_if_newmelted<real_t>()   // functor for r2 of newly melted particles
      );

      // Compute per-cell 3rd moment of newly melted particles (sum of n*r^3). It is stored in count_mom
      moms_all();
      moms_calc(rw2_melted.begin(), count_n, real_t(1.5), true);
      nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd wet moment) of newly melted particles");

      // Copy to melted_mom3 array
      thrust_device::vector<real_t> &melted_mom3(tmp_device_real_cell);
      if (count_n != n_cell)
        thrust::fill(melted_mom3.begin(), melted_mom3.end(), real_t(0.));

      // multiplying by -1 (dth < 0 for melting)
      thrust::transform(
          count_mom.begin(), count_mom.begin() + count_n,
          thrust::make_permutation_iterator(melted_mom3.begin(), count_ijk.begin()),
          thrust::negate<real_t>()
      );

      // update th according to the melted volume per cell
      update_th_freezing(melted_mom3);
    }

  }
}
