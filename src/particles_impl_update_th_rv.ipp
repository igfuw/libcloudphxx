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
      struct dth
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

    // update th and rv according to change in liquid water volume
    // particles have to be sorted
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::update_th_rv(
      thrust_device::vector<real_t> &drv, // change in water volume
      bool update_perSD
    ) 
    {   
      if(!sorted) throw std::runtime_error("update_th_rv called on an unsorted set");

      // multiplying dv*dm_3 by -rho_w*4/3*pi
      thrust::transform(
        drv.begin(), drv.end(),                  // input - 1st arg
        thrust::make_constant_iterator<real_t>(  // input - 2nd arg
          - common::moist_air::rho_w<real_t>() / si::kilograms * si::cubic_metres
          * real_t(4./3) * pi<real_t>()
        ),
        drv.begin(),                             // output
        thrust::multiplies<real_t>()
      );

      // dividing by dv
      thrust::transform(
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()),
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()) + count_n,  // input - 1st arg
        thrust::make_permutation_iterator(dv.begin(),  count_ijk.begin()),            // input - 2nd arg
        thrust::make_permutation_iterator(drv.begin(), count_ijk.begin()),            // output
        thrust::divides<real_t>() 
      ); 

      // dividing d(rhod_rv) by rhod
      thrust::transform(
        thrust::make_permutation_iterator(drv.begin(),  count_ijk.begin()), 
        thrust::make_permutation_iterator(drv.begin(),  count_ijk.begin()) + count_n, // input - 1st arg
        thrust::make_permutation_iterator(rhod.begin(), count_ijk.begin()),           // input - 2nd arg
        thrust::make_permutation_iterator(drv.begin(),  count_ijk.begin()),           // output (in place)
        thrust::divides<real_t>()
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
     
      if(update_perSD)
        thrust::transform(
          sstp_tmp_rv.begin(), sstp_tmp_rv.end(),  // input - 1st arg
          thrust::make_permutation_iterator(drv.begin(), ijk.begin()),           // input - 2nd arg
          sstp_tmp_rv.begin(),            // output
          thrust::plus<real_t>() 
        );

      // updating th
      {
        typedef thrust::zip_iterator<thrust::tuple<
          typename thrust_device::vector<real_t>::iterator,
          typename thrust_device::vector<real_t>::iterator,
          typename thrust_device::vector<real_t>::iterator
        > > zip_it_t;

        thrust_device::vector<real_t> &delta_th(tmp_device_real_cell1);

        // calc dth
        thrust::transform(
          zip_it_t(thrust::make_tuple(  
            drv.begin(),      // 
            T.begin(),        // dth = drv * d_th_d_rv(T, th)
            th.begin()        //
          )),
          zip_it_t(thrust::make_tuple( 
            drv.begin(),      // 
            T.begin(),        // dth = drv * d_th_d_rv(T, th)
            th.begin()        //
          )) + n_cell,
          delta_th.begin(),
          detail::dth<real_t>()     // func
        );

        // apply dth
        thrust::transform(
          th.begin(), th.end(),          // input - 1st arg
          delta_th.begin(),
          th.begin(),                 // output
          thrust::plus<real_t>()
        );
        if(update_perSD)
          thrust::transform(
            sstp_tmp_th.begin(), sstp_tmp_th.end(),          // input - 1st arg
            delta_th.begin(),
            sstp_tmp_th.begin(),                 // output
            thrust::plus<real_t>()
          );
      }
    }
  };  
};
