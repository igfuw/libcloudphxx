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
      template <typename real_t>
      struct dth : thrust::unary_function<const thrust::tuple<real_t, real_t, real_t>&, real_t>
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
    };

    // update th and rv according to change in 3rd specific wet moments
    // particles have to be sorted
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::update_th_rv(
      ref_grid<real_t> &drv // change in water vapor mixing ratio
    ) 
    {   
      if(!sorted) throw std::runtime_error("update_th_rv called on an unsorted set");
      nancheck(drv.get_ref(), "update_th_rv: input drv");

      // multiplying specific 3rd moms diff  by -rho_w*4/3*pi
      thrust::transform(
        drv.begin_ref(), drv.end_ref(),                  // input - 1st arg
        thrust::make_constant_iterator<real_t>(  // input - 2nd arg
          - common::moist_air::rho_w<real_t>() / si::kilograms * si::cubic_metres
          * real_t(4./3) * pi<real_t>()
        ),
        drv.begin_ref(),                             // output
        thrust::multiplies<real_t>()
      );  

      // updating rv 
      assert(*thrust::min_element(rv.begin(), rv.end()) >= 0);
      thrust::transform(
        rv.begin_ref(), rv.end_ref(),  // input - 1st arg
        drv.begin_ref(),           // input - 2nd arg
        rv.begin_ref(),            // output
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
        thrust::transform(
          th.begin(), th.end(),          // input - 1st arg
          thrust::make_transform_iterator(
            zip_it_t(thrust::make_tuple(  
              drv.begin_ref(),      // 
              T.begin_ref(),        // dth = drv * d_th_d_rv(T, th)
              th.begin_ref()        //
            )),
            detail::dth<real_t>()
          ),
          th.begin(),                 // output
          thrust::plus<real_t>()
        );
      }
      nancheck(th, "update_th_rv: th after update");
    }

    // update particle-specific  cell state; done on refined grid
    // particles have to be sorted
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::update_pstate(
      thrust_device::vector<real_t> &pstate, // cell characteristic
      thrust_device::vector<real_t> &pdstate // change in cell characteristic
    ) 
    {   
      if(!sorted) throw std::runtime_error("update_th_rv called on an unsorted set");

      // cell-wise change in state
      ref_grid<real_t> &dstate(tmp_device_real_cell);
      // init dstate with 0s
      thrust::fill(dstate.begin_ref(), dstate.end_ref(), real_t(0));
      // calc sum of pdstate in each cell
      thrust::pair<
        thrust_device::vector<thrust_size_t>::iterator,
        typename thrust_device::vector<real_t>::iterator
      > it_pair = thrust::reduce_by_key(
        sorted_ijk.begin_ref(), sorted_ijk.end_ref(),
        thrust::make_permutation_iterator(pdstate.begin(), sorted_id.begin_ref()),
        count_ijk.begin_ref(),
        count_mom.begin_ref()
      );
      count_n.get_ref() = it_pair.first - count_ijk.begin_ref();

      // add this sum to dstate
      thrust::transform(
        count_mom.begin_ref(), count_mom.begin_ref() + count_n.get_ref(),             // input - 1st arg
        thrust::make_permutation_iterator(dstate.begin_ref(), count_ijk.begin_ref()), // 2nd arg
        thrust::make_permutation_iterator(dstate.begin_ref(), count_ijk.begin_ref()), // output
        thrust::plus<real_t>()
      );

      // add dstate to pstate
      thrust::transform(
        pstate.begin(), pstate.end(),
        thrust::make_permutation_iterator(dstate.begin_ref(), ijk.begin_ref()),
        pstate.begin(),
        thrust::plus<real_t>()
      );
    }

    // update cell state based on particle-specific cell state
    // particles have to be sorted
    /*
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
    */
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::update_state(
      ref_grid<real_t> &state, // cell state
      thrust_device::vector<real_t> &pstate // particle-specific cell state (same for all particles in one cell)
    ) 
    {
     thrust::copy(
       pstate.begin(), pstate.end(),
       thrust::make_permutation_iterator(state.begin_ref(), ijk.begin_ref())
     );   
    }
  };  
};
