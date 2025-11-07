// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <libcloudph++/common/theta_dry.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      // change of th during condensation
      template <typename real_t>
      struct dth //: thrust::unary_function<const thrust::tuple<real_t, real_t, real_t>&, real_t>
      {
        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t, real_t> &tpl) const
        {
          const quantity<si::dimensionless, real_t> 
            drv      = thrust::get<0>(tpl);
          const quantity<si::temperature, real_t> 
            T        = thrust::get<1>(tpl) * si::kelvins;
          const quantity<si::temperature, real_t> 
            th       = thrust::get<2>(tpl) * si::kelvins;

          return drv * common::theta_dry::d_th_d_rv(T, th) / si::kelvins;
        }
      };
      // change of th during sublimation
      template <typename real_t>
      struct dth_subl //: thrust::unary_function<const thrust::tuple<real_t, real_t, real_t>&, real_t>
        {
          BOOST_GPU_ENABLED
          real_t operator()(const thrust::tuple<real_t, real_t, real_t> &tpl) const
          {
            const quantity<si::dimensionless, real_t>
              drv      = thrust::get<0>(tpl);
            const quantity<si::temperature, real_t>
              T        = thrust::get<1>(tpl) * si::kelvins;
            const quantity<si::temperature, real_t>
              th       = thrust::get<2>(tpl) * si::kelvins;

            return drv * common::theta_dry::d_th_d_rv_subl(T, th) / si::kelvins;
          }
        };

      // change of th during freezing
      template <typename real_t>
      struct dth_freezing //: thrust::unary_function<const thrust::tuple<real_t, real_t, real_t>&, real_t>
      {
        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t, real_t> &tpl) const
        {
          const quantity<si::dimensionless, real_t>
            drw      = thrust::get<0>(tpl);
          const quantity<si::temperature, real_t>
            T        = thrust::get<1>(tpl) * si::kelvins;
          const quantity<si::temperature, real_t>
            th       = thrust::get<2>(tpl) * si::kelvins;

          return drw * common::theta_dry::d_th_d_rw_freeze(T, th) / si::kelvins;
        }
      };
    };

    // update th and rv after condensation / sublimation according to change in 3rd specific wet moments
    // particles have to be sorted
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::update_th_rv(
      thrust_device::vector<real_t> &drv, // change in water vapor mixing ratio
      phase_change phase // enum for cond/subl, the default is condensation
    ) 
    {   
      if(!sorted) throw std::runtime_error("libcloudph++: update_th_rv called on an unsorted set");
      nancheck(drv, "update_th_rv: input drv");

      // multiplying specific 3rd moms diff  by -rho_w*4/3*pi
      thrust::transform(
        drv.begin(), drv.end(),                  // input - 1st arg
        thrust::make_constant_iterator<real_t>(  // input - 2nd arg
          - common::moist_air::rho_w<real_t>() / si::kilograms * si::cubic_metres
          * real_t(4./3) * pi<real_t>()
        ),
        drv.begin(),                             // output
        thrust::multiplies<real_t>()
      );  

      // updating rv 
      assert(*thrust::min_element(rv.begin(), rv.end()) >= 0);
      thrust::transform(
        rv.begin(), rv.end(),  // input - 1st arg
        drv.begin(),           // input - 2nd arg
        rv.begin(),            // output
        thrust::plus<real_t>() 
      );
      assert(*thrust::min_element(rv.begin(), rv.end()) >= 0);
      nancheck(rv, "update_th_rv: rv after update");

      // updating th
      {
        typedef thrust::zip_iterator<thrust::tuple<
          typename thrust_device::vector<real_t>::iterator,
          typename thrust_device::vector<real_t>::iterator,
          typename thrust_device::vector<real_t>::iterator
        > > zip_it_t;

        // apply dth
      if (phase == phase_change::condensation)
      {
        thrust::transform(
          th.begin(), th.end(),          // input - 1st arg
          thrust::make_transform_iterator(
            zip_it_t(thrust::make_tuple(
              drv.begin(),      //
              T.begin(),        // dth = drv * d_th_d_rv(T, th)
              th.begin()        //
            )),
            detail::dth<real_t>()
          ),
          th.begin(),                 // output
          thrust::plus<real_t>()
        );
      }
      else if (phase == phase_change::sublimation)
      {
        thrust::transform(
          th.begin(), th.end(),          // input - 1st arg
          thrust::make_transform_iterator(
            zip_it_t(thrust::make_tuple(
              drv.begin(),      //
              T.begin(),        // dth = drv * d_th_d_rv(T, th)
              th.begin()        //
            )),
            detail::dth_subl<real_t>()
          ),
          th.begin(),                 // output
          thrust::plus<real_t>()
        );
        }

      }
      nancheck(th, "update_th_rv: th after update");
    }

    // update th for freezing
    // particles have to be sorted
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::update_th_freezing(
      thrust_device::vector<real_t> &drw // change of specific 3rd moment of liquid per cell
    )
    {
      if(!sorted) throw std::runtime_error("libcloudph++: update_th_freezing called on an unsorted set");
      nancheck(drw, "update_th_freezing: input drw");

      // Calculating the change of liquid mixing ratio per cell (multiplying specific 3rd mom by rho_w*4/3*pi)
      thrust::transform(
        drw.begin(), drw.end(),   // input - 1st arg
        thrust::make_constant_iterator<real_t>(  // input - 2nd arg
          common::moist_air::rho_w<real_t>() / si::kilograms * si::cubic_metres
          * real_t(4./3) * pi<real_t>()
        ),
        drw.begin(),                             // output
        thrust::multiplies<real_t>()
      );

      // updating th
      {
        typedef thrust::zip_iterator<thrust::tuple<
          typename thrust_device::vector<real_t>::iterator,
          typename thrust_device::vector<real_t>::iterator,
          typename thrust_device::vector<real_t>::iterator
        > > zip_it_t;

        // apply dth
        thrust::transform(
          th.begin(), th.end(),          // input - 1st arg
          thrust::make_transform_iterator(
            zip_it_t(thrust::make_tuple(
              drw.begin(),
              T.begin(),
              th.begin()
            )),
            detail::dth_freezing<real_t>()
          ),
          th.begin(),                 // output
          thrust::plus<real_t>()
        );
      }
      nancheck(th, "update_th_freezing: th after update");
    }

    // update particle-specific cell state
    // particles have to be sorted
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::update_pstate(
      thrust_device::vector<real_t> &pstate, // cell characteristic
      thrust_device::vector<real_t> &pdstate // change in cell characteristic
    ) 
    {   
      if(!sorted) throw std::runtime_error("libcloudph++: update_uh_rv called on an unsorted set");

      // cell-wise change in state
      auto dstate_g = tmp_device_real_cell.get_guard();
      thrust_device::vector<real_t> &dstate(dstate_g.get());
      // init dstate with 0s
      thrust::fill(dstate.begin(), dstate.end(), real_t(0));
      // calc sum of pdstate in each cell
      thrust::pair<
        thrust_device::vector<thrust_size_t>::iterator,
        typename thrust_device::vector<real_t>::iterator
      > it_pair = thrust::reduce_by_key(
        sorted_ijk.begin(), sorted_ijk.end(),
        thrust::make_permutation_iterator(pdstate.begin(), sorted_id.begin()),
        count_ijk.begin(),
        count_mom.begin()
      );
      count_n = it_pair.first - count_ijk.begin();

      // add this sum to dstate
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,                    // input - 1st arg
        thrust::make_permutation_iterator(dstate.begin(), count_ijk.begin()), // 2nd arg
        thrust::make_permutation_iterator(dstate.begin(), count_ijk.begin()), // output
        thrust::plus<real_t>()
      );

      // add dstate to pstate
      thrust::transform(
        pstate.begin(), pstate.end(),
        thrust::make_permutation_iterator(dstate.begin(), ijk.begin()),
        pstate.begin(),
        thrust::plus<real_t>()
      );
    }

    // update cell state based on particle-specific cell state
    // particles have to be sorted
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::update_state(
      thrust_device::vector<real_t> &state, // cell state
      thrust_device::vector<real_t> &pstate // particle-specific cell state (same for all particles in one cell)
    ) 
    {
     thrust::copy(
       pstate.begin(), pstate.end(),
       thrust::make_permutation_iterator(state.begin(), ijk.begin())
     );   
    }
  };  
};
