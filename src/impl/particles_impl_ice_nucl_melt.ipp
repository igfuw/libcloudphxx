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
      // Singular immersion freezing (Shima et al., 2020)
      // Functor to update ice flag, rw2, a, c, rho_i, of frozen droplets
      template<class real_t>
      class singular_freeze
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

          if (T_freeze >=  T && RH >= real_t(1)) // condition for freezing
          {
            ice = real_t(1);
            rw2  = real_t(0);
            rho_i = common::moist_air::rho_i<real_t>().value();
            a   = pow(rw2, real_t(0.5)) * pow(common::moist_air::rho_w<real_t>() / common::moist_air::rho_i<real_t>(), real_t(1./3.));
            c   = pow(rw2, real_t(0.5)) * pow(common::moist_air::rho_w<real_t>() / common::moist_air::rho_i<real_t>(), real_t(1./3.));
          }
        }
      };

      // Time-dependent immersion freezing (Arabas et al., 2025)
      // Functor to update ice flag, rw2, a, c, rho_i, of frozen droplets
      template<class real_t>
      class time_dep_freeze
      {
        const real_t dt;
        //const common::ice_nucleation::INP_t INP_type;

      public:
        BOOST_GPU_ENABLED
        time_dep_freeze(const real_t &dt) : dt(dt)
        {}

        BOOST_GPU_ENABLED
        void operator()(thrust::tuple<
            real_t&, real_t&, real_t&, real_t&, real_t&, // to be updated (ice, rw2, a, c, rho_i)
            const real_t&, const real_t&, const real_t& // rd2_insol, u01, T
          > tpl) const
        {
          auto& ice   = thrust::get<0>(tpl);
          auto& rw2   = thrust::get<1>(tpl);
          auto& a     = thrust::get<2>(tpl);
          auto& c     = thrust::get<3>(tpl);
          auto& rho_i = thrust::get<4>(tpl);

          const real_t rd2_insol = thrust::get<5>(tpl);
          const real_t u01 = thrust::get<6>(tpl);
          const real_t T  = thrust::get<7>(tpl);

          if (u01 < common::ice_nucleation::p_freeze<real_t>(common::ice_nucleation::INP_t::mineral, rd2_insol, T, dt))
          {
            ice = real_t(1);
            rw2  = real_t(0);
            rho_i = common::moist_air::rho_i<real_t>().value();
            a   = pow(rw2, real_t(0.5)) * pow(common::moist_air::rho_w<real_t>() / common::moist_air::rho_i<real_t>(), real_t(1./3.));
            c   = pow(rw2, real_t(0.5)) * pow(common::moist_air::rho_w<real_t>() / common::moist_air::rho_i<real_t>(), real_t(1./3.));
          }
        }
      };

      // Functor to update ice flag, rw2, a, c, rho_i, of melted ice
      template<class real_t>
      class melt
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

          if (thrust::get<5>(tpl) > real_t(273.15)) // if T > 0 C
          {
            ice = real_t(0);
            rw2  = pow(common::moist_air::rho_i<real_t>() / common::moist_air::rho_w<real_t>() * c , real_t(2./3.)) * pow(a , real_t(4./3.));
            rho_i = real_t(0);
            a   = real_t(0);
            c   = real_t(0);
          }
        }
      };

    }

    // Immersion freezing and melting
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::ice_nucl_melt(const real_t &dt, const bool time_dep_ice_nucl) {

      hskpng_sort();

      // A vector to store liquid 3rd moments
      thrust_device::vector<real_t> &drw(tmp_device_real_cell);

      // Compute per-cell 3rd moment of liquid droplets (sum of n*r^3) before freezing/melting. It is stored in count_mom
      moms_eq0(ice.begin()); // choose particles with ice=0
      moms_calc(rw2.begin(), real_t(1.5));
      nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd wet moment) of droplets before freezing/melting");
      if (count_n != n_cell)
        thrust::fill(drw.begin(), drw.end(), real_t(0.));

      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,
        thrust::make_permutation_iterator(drw.begin(), count_ijk.begin()),
        thrust::negate<real_t>()
      );

      // Change liquid droplets to ice under the freezing condition

      if (time_dep_ice_nucl) // time dependent freezing based on Arabas et al., 2025
      {
        rand_u01(n_part); // random numbers between [0,1] for each particle
        thrust::for_each(
          thrust::make_zip_iterator(thrust::make_tuple(
            ice.begin(),
            rw2.begin(),
            ice_a.begin(),
            ice_c.begin(),
            ice_rho.begin(),
            rd2_insol.begin(),
            u01.begin(),
            thrust::make_permutation_iterator(T.begin(), ijk.begin())
          )),
          thrust::make_zip_iterator(thrust::make_tuple(
            ice.begin(),
            rw2.begin(),
            ice_a.begin(),
            ice_c.begin(),
            ice_rho.begin(),
            rd2_insol.begin(),
            u01.begin(),
            thrust::make_permutation_iterator(T.begin(), ijk.begin())
          )) + n_part,
            detail::time_dep_freeze<real_t>(dt)  // functor for updating (ice, rw2, a, c, rho_i) if freezing condition satisfied
        );
      }
      else  // singular freezing based on Shima et al., 2020
      {
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
            detail::singular_freeze<real_t>()  // functor for updating (ice, rw2, a, c, rho_i) if freezing condition satisfied
        );
      }

      // Change ice to liquid droplets under the melting condition
      thrust::for_each(
        thrust::make_zip_iterator(thrust::make_tuple(
          ice.begin(),
          rw2.begin(),
          ice_a.begin(),
          ice_c.begin(),
          ice_rho.begin(),
          thrust::make_permutation_iterator(T.begin(), ijk.begin())
        )),
        thrust::make_zip_iterator(thrust::make_tuple(
          ice.begin(),
          rw2.begin(),
          ice_a.begin(),
          ice_c.begin(),
          ice_rho.begin(),
          thrust::make_permutation_iterator(T.begin(), ijk.begin())
        )) + n_part,
          detail::melt<real_t>()  // functor for updating (ice, rw2, a, c, rho_i) if melting condition satisfied
      );

      // Compute per-cell 3rd moment of liquid droplets after freezing/melting. It is stored in count_mom
      moms_eq0(ice.begin()); // choose particles with ice=0
      moms_calc(rw2.begin(), real_t(1.5));
      nancheck_range(count_mom.begin(), count_mom.begin() + count_n, "count_mom (3rd wet moment) of droplets after freezing/melting");

      // Compute the difference between liquid moments before and after
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
