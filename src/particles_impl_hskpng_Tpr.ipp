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
        __device__ real_t operator()(const real_t &th, const real_t &rhod)
       {   
         return common::theta_dry::T<real_t>(
           th   * si::kelvins,
           rhod * si::kilograms / si::cubic_metres
         ) / si::kelvins;
       }   
      }; 

      template <typename real_t>
      struct common__theta_dry__p 
      {
        __device__ real_t operator()(const thrust::tuple<real_t, real_t, real_t> &tpl)
       {   
         const real_t rhod = thrust::get<0>(tpl);
         const real_t rv = thrust::get<1>(tpl);
         const real_t T = thrust::get<2>(tpl);

         return common::theta_dry::p<real_t>(
           rhod * si::kilograms / si::cubic_metres,
           rv,
           T * si::kelvins 
         ) / si::pascals;
       }   
      }; 

      template <typename real_t>
      struct RH
      {   
        __device__ 
        real_t operator()(const thrust::tuple<real_t, real_t, real_t> &tpl) 
        {
          const real_t rhod = thrust::get<0>(tpl);
          const real_t rv = thrust::get<1>(tpl);
          const real_t T = thrust::get<2>(tpl);

          return 
	    (rhod * rv * si::kilograms / si::cubic_metres)
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

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_Tpr()
    {   
      using namespace thrust::placeholders;

      // T  = common::theta_dry::T<real_t>(th, rhod);
      thrust::transform(
        th.begin(), th.end(),      // input - first arg
        rhod.begin(),              // input - second arg
        T.begin(),                 // output
        detail::common__theta_dry__T<real_t>() 
      );

      {
        typedef thrust::zip_iterator<
          thrust::tuple<
            typename thrust_device::vector<real_t>::iterator,
            typename thrust_device::vector<real_t>::iterator,
            typename thrust_device::vector<real_t>::iterator
          >
        > zip_it_t;

        // p  = common::theta_dry::p<real_t>(rhod, r, T); 
        thrust::transform(
          zip_it_t(thrust::make_tuple(rhod.begin(), rv.begin(), T.begin())), // input - begin
          zip_it_t(thrust::make_tuple(rhod.end(),   rv.end(),   T.end()  )), // input - end
          p.begin(),                                                         // output
          detail::common__theta_dry__p<real_t>()
        );

        // RH = p_v / p_vs = rhod * rv * R_v * T / p_vs
	thrust::transform(
	  zip_it_t(thrust::make_tuple(rhod.begin(), rv.begin(), T.begin())),  // input - begin
	  zip_it_t(thrust::make_tuple(rhod.end(),   rv.end(),   T.end()  )),  // input - end
	  RH.begin(),                                                         // output
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
