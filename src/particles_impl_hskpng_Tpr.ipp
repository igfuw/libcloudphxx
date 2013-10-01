// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <libcloudph++/common/theta_dry.hpp>
#include <libcloudph++/common/vterm.hpp> // TODO: should be viscosity!

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t>
      struct common__theta_dry__T 
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
      struct common__theta_dry__p 
      {
        __device__ real_t operator()(const real_t &rhod_T, const real_t &r)
       {   
         return common::theta_dry::p<real_t>(
           rhod_T * si::kilograms / si::cubic_metres * si::kelvins,
           r 
         ) / si::pascals;
       }   
      }; 

      template <typename real_t>
      struct RH
      {   
        __device__ 
        real_t operator()(const thrust::tuple<real_t, real_t> &tpl) 
        {
          const real_t rhod_rv = thrust::get<0>(tpl);
          const real_t T = thrust::get<1>(tpl);

          return 
	    (rhod_rv * si::kilograms / si::cubic_metres)
	    * common::moist_air::R_v<real_t>()
	    * (T * si::kelvins)
	    / common::const_cp::p_vs(T * si::kelvins);
        }
      }; 

      
      template <typename real_t>
      struct common__vterm__visc // TODO: rename it! (vterm) visc_eta?
      {
        __device__
        real_t operator()(const real_t &T)
        {
          return common::vterm::visc(T * si::kelvins) / si::pascals / si::seconds;
        }
      };
    };

    template <typename real_t, int device>
    void particles_t<real_t, device>::impl::hskpng_Tpr()
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
        // TODO: rewrite with zip iterator?
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

      // RH = p_v / p_vs = rhod_rv * R_v * T / p_vs
      {
        typedef thrust::zip_iterator<
          thrust::tuple<
            typename thrust_device::vector<real_t>::iterator,
            typename thrust_device::vector<real_t>::iterator
          >
        > zip_it_t;

	thrust::transform(
	  zip_it_t(thrust::make_tuple(rhod_rv.begin(), T.begin())),  // input - begin
	  zip_it_t(thrust::make_tuple(rhod_rv.end(),   T.end()  )),  // input - end
	  RH.begin(),                                                // output
	  detail::RH<real_t>()
	);
      }
 
      // dynamic viscosity
      {
        thrust::transform(
          T.begin(), T.end(), // 1st arg
          eta.begin(),        // output
          detail::common__vterm__visc<real_t>()
        );
      }
    }
  };  
};
