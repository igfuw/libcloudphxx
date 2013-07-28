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
      struct common__theta_dry__T /* : thrust::binary_function<real_t, real_t, real_t> */ 
      {
        __device__ real_t operator()(const real_t &rhod_th, const real_t &rhod)
       {   
         return common::theta_dry::T<real_t>(
           rhod_th * si::kilograms / si::cubic_metres * si::kelvins,
           rhod    * si::kilograms / si::cubic_metres
         ) / si::kelvins;
       }   
      }; 

      template <typename real_t>
      struct common__theta_dry__p /* : thrust::binary_function<real_t, real_t, real_t> */ 
      {
        __device__ real_t operator()(const real_t &rhod_T, const real_t &r)
       {   
         return common::theta_dry::p<real_t>(
           rhod_T * si::kilograms / si::cubic_metres * si::kelvins,
           r 
         ) / si::pascals;
       }   
      }; 
    };

    template <typename real_t, int device>
    void particles<real_t, device>::impl::hskpng_Tpr()
    {   
      using namespace thrust::placeholders;

      // r  = rhod_rv / rhod;
      thrust::transform(
	rhod_rv.begin(), rhod_rv.end(), // input - first arg
	rhod.begin(),                   // input - second arg
	r.begin(),                      // output
	_1 / _2 
      );

      // T  = common::theta_dry::T<real_t>(rhod_th, rhod);
      thrust::transform(
        rhod_th.begin(), rhod_th.end(), // input - first arg
        rhod.begin(),                   // input - second arg
        T.begin(),                      // output
        detail::common__theta_dry__T<real_t>() 
      );

      // p  = common::theta_dry::p<real_t>(rhod, r, T); 
      {
        thrust_device::vector<real_t> &rhod_T(p); 
        thrust::transform(
          rhod.begin(), rhod.end(),     // input - first arg
          T.begin(),                    // input - second arg
          rhod_T.begin(),               // output
          _1 * _2
        );
        thrust::transform(
          rhod_T.begin(), rhod_T.end(), // input - first arg
          r.begin(),                    // input - second arg
          p.begin(),                    // output (here making it in-place as rhod_T points to p)
          detail::common__theta_dry__p<real_t>()
        );
      }
    }
  };  
};
