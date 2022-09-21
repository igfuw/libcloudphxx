// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#if defined(__NVCC__)
  #include <nvfunctional>
#endif
#include <libcloudph++/common/theta_dry.hpp>
#include <libcloudph++/common/vterm.hpp> // TODO: should be viscosity!
#include <libcloudph++/common/tetens.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
#if defined(__NVCC__)
      using nvstd::function;
#else
      using std::function;
#endif
      template <typename real_t>
      struct common__theta_dry__T_rhod 
      {
       BOOST_GPU_ENABLED 
       real_t operator()(const real_t &th, const real_t &rhod)
       {   
         return common::theta_dry::T<real_t>(
           th   * si::kelvins,
           rhod * si::kilograms / si::cubic_metres
         ) / si::kelvins;
       }   
      }; 

      template <typename real_t>
      struct common__theta_dry__T_p 
      {
       template <class tpl_t>
       BOOST_GPU_ENABLED 
       real_t operator()(const real_t &th, const tpl_t &tpl) // tpl: (rv, p)
       {   
         return th * 
                common::theta_std::exner<real_t>(thrust::get<1>(tpl)  * si::pascals);
       }   
      }; 

      template <typename real_t>
      struct common__theta_dry__p : thrust::unary_function<const thrust::tuple<real_t, real_t, real_t>&, real_t>
      {
       BOOST_GPU_ENABLED
       real_t operator()(const thrust::tuple<real_t, real_t, real_t> &tpl)
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
      BOOST_GPU_ENABLED
      real_t RH_pv_cc(const real_t &p, const real_t &rv, const real_t &T)
      {
        return real_t( 
          common::moist_air::p_v(p * si::pascals, quantity<si::dimensionless, real_t>(rv))
          / common::const_cp::p_vs(T * si::kelvins));
      }

      template <typename real_t>
      BOOST_GPU_ENABLED
      real_t RH_rv_cc(const real_t &p, const real_t &rv, const real_t &T)
      {
        return real_t(
          rv
          / common::const_cp::r_vs(T * si::kelvins, p * si::pascals)); 
      }

      template <typename real_t>
      BOOST_GPU_ENABLED
      real_t RH_pv_tet(const real_t &p, const real_t &rv, const real_t &T)
      {
        return real_t(
          common::moist_air::p_v(p * si::pascals, quantity<si::dimensionless, real_t>(rv))
          / common::tetens::p_vs(T * si::kelvins)); 
      }

      template <typename real_t>
      BOOST_GPU_ENABLED
      real_t RH_rv_tet(const real_t &p, const real_t &rv, const real_t &T)
      {
        return real_t( 
          rv
          / common::tetens::r_vs(T * si::kelvins, p * si::pascals)); 
      }

      template <typename real_t>
      struct RH : thrust::unary_function<const thrust::tuple<real_t, real_t, real_t>&, real_t>
      {   
        /*
 *      on CUDA 8.0 the code below compiles, but gives memory errors at runtime
 *      probably due to some __host__ / __device__ code mismatch
 *      TODO: fix it
 */

/* 
        const function<real_t(const real_t&, const real_t&, const real_t&)> RH_fun;

        // the type of formula to be used for RH
        RH(RH_formula_t::RH_formula_t RH_formula):
          RH_fun(
            RH_formula == RH_formula_t::pv_cc  ? RH_pv_cc<real_t>  : 
            RH_formula == RH_formula_t::rv_cc  ? RH_rv_cc<real_t>  : 
            RH_formula == RH_formula_t::pv_tet ? RH_pv_tet<real_t> :
            RH_formula == RH_formula_t::rv_tet ? RH_rv_tet<real_t> :
            RH_pv_cc<real_t> // if unrecognised option, use pv_cc, TODO: throw assert?
          )
        {}

        BOOST_GPU_ENABLED 
        real_t operator()(const thrust::tuple<real_t, real_t, real_t> &tpl) 
        {
          return RH_fun(thrust::get<0>(tpl), thrust::get<1>(tpl), thrust::get<2>(tpl)); // p, rv, T
        }
*/

        // an alternative implementation with formula choice at functor call
        const RH_formula_t RH_formula;
        // the type of formula to be used for RH
        RH(RH_formula_t RH_formula):
          RH_formula(RH_formula)
        {}

        BOOST_GPU_ENABLED 
        real_t operator()(const thrust::tuple<real_t, real_t, real_t> &tpl)  // p, rv, T
        {
          switch (RH_formula)
          {
            case RH_formula_t::pv_cc:
              return RH_pv_cc<real_t>(thrust::get<0>(tpl), thrust::get<1>(tpl), thrust::get<2>(tpl));
            case RH_formula_t::rv_cc:
              return RH_rv_cc<real_t>(thrust::get<0>(tpl), thrust::get<1>(tpl), thrust::get<2>(tpl));
            case RH_formula_t::pv_tet:
              return RH_pv_tet<real_t>(thrust::get<0>(tpl), thrust::get<1>(tpl), thrust::get<2>(tpl));
            case RH_formula_t::rv_tet:
              return RH_rv_tet<real_t>(thrust::get<0>(tpl), thrust::get<1>(tpl), thrust::get<2>(tpl));
            default:
              assert(0);
              return 0.;
          }
        }
      }; 
      
      template <typename real_t>
      struct common__vterm__visc : thrust::unary_function<const real_t&, real_t>// TODO: rename it! (vterm) visc_eta?
      {
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &T)
        {
          return common::vterm::visc(T * si::kelvins) / si::pascals / si::seconds;
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_Tpr()
    {   
      if(!const_p) // variable pressure
      {
        // T  = common::theta_dry::T<real_t>(th, rhod);
        thrust::transform(
          th.begin_ref(), th.end_ref(),  // input - first arg
          thrust::make_permutation_iterator(rhod.begin(), ijk_ref2ijk.begin()),             // input - second arg
          T.begin_ref(),                 // output
          detail::common__theta_dry__T_rhod<real_t>() 
        );
      }
      else // external pressure profile
      {
        // T = th * exner(p_tot)
        thrust::transform(
          th.begin_ref(), th.end_ref(),      // input - first arg
          thrust::make_zip_iterator(thrust::make_tuple(rv.begin_ref(), p.begin_ref())), // input - second and third args
          T.begin_ref(),                 // output
          detail::common__theta_dry__T_p<real_t>() 
        );
      }

      {
        typedef thrust::zip_iterator<
          thrust::tuple<
            typename thrust_device::vector<real_t>::iterator,
            typename thrust_device::vector<real_t>::iterator,
            typename thrust_device::vector<real_t>::iterator
          >
        > zip_it_t;

        if(!const_p)
        {
          // p  = common::theta_dry::p<real_t>(rhod, r, T); 
          zip_it_t it(thrust::make_tuple(
            thrust::make_permutation_iterator(rhod.begin(), ijk_ref2ijk.begin()), 
            rv.begin_ref(), 
            T.begin_ref()
          ));
          thrust::transform(
            it,                    // input - begin
            it + n_cell.get_ref(), // input - end
            p.begin_ref(),         // output
            detail::common__theta_dry__p<real_t>()
          );
        }

        thrust::transform(
          zip_it_t(thrust::make_tuple(p.begin_ref(), rv.begin_ref(), T.begin_ref())),  // input - begin
          zip_it_t(thrust::make_tuple(p.end_ref(),   rv.end_ref(),   T.end_ref()  )),  // input - end
          RH.begin_ref(),                                                  // output
          detail::RH<real_t>(opts_init.RH_formula)
        );
      }
 
      // dynamic viscosity
      {
        // on refined grid
        thrust::transform(
          T.begin_ref(), T.end_ref(), // 1st arg
          eta.begin_ref(),        // output
          detail::common__vterm__visc<real_t>()
        );
        // on normal grid: how if we dont know T, because we dont know th?
        // TODO: require input of normal th? 
        // or use averaging from refined to normal?
        // cumbersome: 
        // 0. make copies of eta and ijk_ref2ijk
        // 1. sort copies of eta by ijk_ref2ijk
        // 2. reduce_by_key and add the result to normal eta
      }


      // adjusting dv if using a parcel set-up (1kg of dry air)
      if (n_dims == 0)
      {
        namespace arg = thrust::placeholders;
        thrust::transform(
          rhod.begin(), rhod.end(), // input
          dv.begin(),               // output
          real_t(1) / arg::_1
        );
      }
    }
  };  
};
