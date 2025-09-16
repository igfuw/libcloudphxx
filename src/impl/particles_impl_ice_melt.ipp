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

      // Functor to update ice flag, rw2, a, c, rho_i, of melted ice
      template<class real_t>
      class melting_update
      {
      public:
        BOOST_GPU_ENABLED
        void operator()(thrust::tuple<
            real_t&, real_t&, real_t&, real_t&, real_t&, // to be updated (ice, rw2, a, c, rho_i)
            const real_t& // ambient T
          > tpl) const
        {
          auto& ice   = thrust::get<0>(tpl);
          auto& rw2   = thrust::get<1>(tpl);
          auto& a     = thrust::get<2>(tpl);
          auto& c     = thrust::get<3>(tpl);
          auto& rho_i = thrust::get<4>(tpl);

          if (detail::melting_cond<real_t>()(thrust::get<5>(tpl)))
          {
            ice = real_t(0);
            rw2  = pow(common::moist_air::rho_i<real_t>() / common::moist_air::rho_w<real_t>() * c , real_t(2./3.)) * pow(a , real_t(4./3.));
            rho_i = real_t(0);
            a   = real_t(0);
            c   = real_t(0);
          }
        }
      };

    };

    // Melting
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::ice_melt() {

     hskpng_sort();

      // A vector to store liquid 3rd moments
      thrust_device::vector<real_t> &drw(tmp_device_real_cell);
      if (count_n != n_cell)
        thrust::fill(drw.begin(), drw.end(), real_t(0.));

      // Compute per-cell 3rd moment of liquid droplets (sum of n*r^3) before melting. It is stored in count_mom
      moms_eq0(ice.begin()); // choose particles with ice=0
      moms_calc(rw2.begin(), real_t(1.5));
      nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd wet moment) of droplets before melting");

      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,
        thrust::make_permutation_iterator(drw.begin(), count_ijk.begin()),
        thrust::negate<real_t>()
      );

      // Change ice to liquid droplets under the melting condition
      thrust::for_each(
        thrust::make_zip_iterator(thrust::make_tuple(
          ice.begin(),
          rw2.begin(),
          a_ice.begin(),
          c_ice.begin(),
          rho_i.begin(),
          thrust::make_permutation_iterator(T.begin(), ijk.begin())
        )),
        thrust::make_zip_iterator(thrust::make_tuple(
          ice.begin(),
          rw2.begin(),
          a_ice.begin(),
          c_ice.begin(),
          rho_i.begin(),
          thrust::make_permutation_iterator(T.begin(), ijk.begin())
        )) + n_part,
          detail::melting_update<real_t>()  // functor for updating (ice, rw2, a, c, rho_i) if melting condition satisfied
      );

      // Compute per-cell 3rd moment of liquid droplets after melting. It is stored in count_mom
      moms_eq0(ice.begin()); // choose particles with ice=0
      moms_calc(rw2.begin(), real_t(1.5));
      nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd wet moment) of droplets after melting");

      // Compute the difference between liquid moments before and after melting
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,
        thrust::make_permutation_iterator(drw.begin(), count_ijk.begin()),
        thrust::make_permutation_iterator(drw.begin(), count_ijk.begin()),
        thrust::plus<real_t>()
      );

      // Update th according to the melted volume per cell
        update_th_freezing(drw);

    }

  }
}
