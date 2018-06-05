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
         return common::theta_dry::dry2std(th * si::kelvins, quantity<si::dimensionless, real_t>(thrust::get<0>(tpl))) / si::kelvins * 
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
      struct RH : thrust::unary_function<const thrust::tuple<real_t, real_t, real_t>&, real_t>
      {   
#if defined(__NVCC__)
        using nvstd::function
#else
        using std::function
#endif
        const function<real_t(const real_t&, const real_t&, const real_t&)> RH_fun;

        // the type of formula to be used for RH
        RH(RH_formula_t::RH_formula_t RH_formula):
          RH_fun(
            RH_formula == RH_formula_t::pv_cc ? 
              [](const real_t &p, const real_t &rv, const real_t &T){
                return real_t( 
                  common::moist_air::p_v(p * si::pascals, quantity<si::dimensionless, real_t>(rv))
                  / common::const_cp::p_vs(T * si::kelvins));} :

            RH_formula == RH_formula_t::rv_cc ? 
              [](const real_t &p, const real_t &rv, const real_t &T){
                return real_t(
                  rv
                  / common::const_cp::r_vs(T * si::kelvins, p * si::pascals));} :

            RH_formula == RH_formula_t::pv_tet ? 
              [](const real_t &p, const real_t &rv, const real_t &T){
                return real_t(
                  common::moist_air::p_v(p * si::pascals, quantity<si::dimensionless, real_t>(rv))
                  / common::tetens::p_vs(T * si::kelvins));} :

            RH_formula == RH_formula_t::rv_tet ? 
              [](const real_t &p, const real_t &rv, const real_t &T){
                return real_t( 
                  rv
                  / common::tetens::r_vs(T * si::kelvins, p * si::pascals));} :

              [](const real_t &p, const real_t &rv, const real_t &T){
                assert(0);
                return real_t(0.);}
          )
        {}

        BOOST_GPU_ENABLED 
        real_t operator()(const thrust::tuple<real_t, real_t, real_t> &tpl) 
        {
          return RH_fun(thrust::get<0>(tpl), thrust::get<1>(tpl), thrust::get<2>(tpl)); // p, rv, T
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
          th.begin(), th.end(),      // input - first arg
          rhod.begin(),              // input - second arg
          T.begin(),                 // output
          detail::common__theta_dry__T_rhod<real_t>() 
        );
      }
      else // external pressure profile
      {
        // T = dry2std(th_d, rv) * exner(p_tot)
        thrust::transform(
          th.begin(), th.end(),      // input - first arg
          thrust::make_zip_iterator(thrust::make_tuple(rv.begin(), p.begin())), // input - second and third args
          T.begin(),                 // output
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
          thrust::transform(
            zip_it_t(thrust::make_tuple(rhod.begin(), rv.begin(), T.begin())), // input - begin
            zip_it_t(thrust::make_tuple(rhod.end(),   rv.end(),   T.end()  )), // input - end
            p.begin(),                                                         // output
            detail::common__theta_dry__p<real_t>()
          );
        }

        thrust::transform(
          zip_it_t(thrust::make_tuple(p.begin(), rv.begin(), T.begin())),  // input - begin
          zip_it_t(thrust::make_tuple(p.end(),   rv.end(),   T.end()  )),  // input - end
          RH.begin(),                                                      // output
          detail::RH<real_t>(opts_init.RH_formula)
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
