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

        // Functor to update ice flag, rw2, a, c, rho_i, of frozen droplets
        template<class real_t>
        class freezing_update
        {
        public:
          BOOST_GPU_ENABLED
          void operator()(thrust::tuple<
              real_t&, real_t&, real_t&, real_t&, real_t&, // to be updated (ice, rw2, a, c, rho_i)
              const real_t&, const real_t&, const real_t& // T_freeze, T, RH
            > tpl) const
          {
            auto& ice   = thrust::get<0>(tpl);
            auto& rw2   = thrust::get<1>(tpl);
            auto& a     = thrust::get<2>(tpl);
            auto& c     = thrust::get<3>(tpl);
            auto& rho_i = thrust::get<4>(tpl);

            const real_t T_freeze = thrust::get<5>(tpl);
            const real_t T        = thrust::get<6>(tpl);
            const real_t RH       = thrust::get<7>(tpl);

            if (detail::immersion_freeze_cond<real_t>()(thrust::make_tuple(T_freeze, T, RH)))
            {
              ice = real_t(1);
              rw2  = real_t(0);
              rho_i = common::moist_air::rho_i<real_t>().value();
              a   = pow(rw2, real_t(0.5)) * pow(common::moist_air::rho_w<real_t>() / common::moist_air::rho_i<real_t>(), real_t(1./3.));
              c   = pow(rw2, real_t(0.5)) * pow(common::moist_air::rho_w<real_t>() / common::moist_air::rho_i<real_t>(), real_t(1./3.));
            }
          }
        };

         }

    // Immersion freezing
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::ice_nucl() {

      hskpng_sort();

      // A vector to store liquid 3rd moments
      thrust_device::vector<real_t> &drw(tmp_device_real_cell);
      if (count_n != n_cell)
        thrust::fill(drw.begin(), drw.end(), real_t(0.));

      // Compute per-cell 3rd moment of liquid droplets (sum of n*r^3) before freezing. It is stored in count_mom
      moms_eq0(ice.begin()); // choose particles with ice=0
      moms_calc(rw2.begin(), real_t(1.5));
      nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd wet moment) of droplets before freezing");

      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,
        thrust::make_permutation_iterator(drw.begin(), count_ijk.begin()),
        thrust::negate<real_t>()
      );

      // Change liquid droplets to ice under the freezing condition
      thrust::for_each(
        thrust::make_zip_iterator(thrust::make_tuple(
          ice.begin(),
          rw2.begin(),
          ice_a.begin(),
          ice_c.begin(),
          ice_rho.begin(),
          T_freeze.begin(),
          thrust::make_permutation_iterator(T.begin(), ijk.begin()),
          thrust::make_permutation_iterator(RH.begin(), ijk.begin())
        )),
        thrust::make_zip_iterator(thrust::make_tuple(
          ice.begin(),
          rw2.begin(),
          ice_a.begin(),
          ice_c.begin(),
          ice_rho.begin(),
          T_freeze.begin(),
          thrust::make_permutation_iterator(T.begin(), ijk.begin()),
          thrust::make_permutation_iterator(RH.begin(), ijk.begin())
        )) + n_part,
          detail::freezing_update<real_t>()  // functor for updating (ice, rw2, a, c, rho_i) if freezing condition satisfied
      );

      // Compute per-cell 3rd moment of liquid droplets after freezing. It is stored in count_mom
      moms_eq0(ice.begin()); // choose particles with ice=0
      moms_calc(rw2.begin(), real_t(1.5));
      nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd wet moment) of droplets after freezing");

      // Compute the difference between liquid moments before and after freezing
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,
        thrust::make_permutation_iterator(drw.begin(), count_ijk.begin()),
        thrust::make_permutation_iterator(drw.begin(), count_ijk.begin()),
        thrust::plus<real_t>()
      );

      // Update th according to the frozen volume per cell
        update_th_freezing(drw);

    }

  }
}
