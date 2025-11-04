// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <libcloudph++/common/vterm.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t>
      struct get_vt0_bin
      {
        const real_t dlnr, ln_r_min, ln_r_max;
        const int n_bin;

        get_vt0_bin(const real_t &ln_r_min, const real_t &ln_r_max, const int &n_bin):
          ln_r_min(ln_r_min),
          ln_r_max(ln_r_max),
          n_bin(n_bin),
          dlnr((ln_r_max - ln_r_min) / n_bin) {}

        BOOST_GPU_ENABLED
        int operator()(const real_t &rw2)
        {
          real_t lnr = .5 * log(rw2);
          return lnr <= ln_r_min ? 0 :
                 lnr >= ln_r_max ? n_bin-1:
                 (lnr - ln_r_min) / dlnr;
        }
      };
      template <typename real_t>
      struct common__vterm__vt
      {
        vt_t vt_eq; //type of terminal velocity formula to use 

        //ctor
        common__vterm__vt(const vt_t &vt_eq): vt_eq(vt_eq) {}

        BOOST_GPU_ENABLED 
        real_t operator()(
          const real_t &rw2, 
          const thrust::tuple<real_t, real_t, real_t, real_t> &tpl // (T, p, rhod, eta)
        ) {   
#if !defined(__NVCC__)
          using std::sqrt;
#endif

          real_t vt;

          // TODO: move the formula selection to the functor's constructor to avoid the switch in each function call,
          //       see how it is done in RH calculation in hskpng_Tpr
          switch(vt_eq)
          {
            case(vt_t::beard76):
              vt = common::vterm::vt_beard76(
                sqrt(rw2)           * si::metres, // TODO: consider caching rw?
                thrust::get<0>(tpl) * si::kelvins,
                thrust::get<1>(tpl) * si::pascals,
                thrust::get<2>(tpl) * si::kilograms / si::cubic_metres,
                thrust::get<3>(tpl) * si::pascals * si::seconds
              ) / si::metres_per_second;
              break;

            case(vt_t::beard77):
              vt = 
                common::vterm::vt_beard77_fact(
                  sqrt(rw2)           * si::metres, // TODO: consider caching rw?
                  thrust::get<1>(tpl) * si::pascals,
                  thrust::get<2>(tpl) * si::kilograms / si::cubic_metres,
                  thrust::get<3>(tpl) * si::pascals * si::seconds
                ) * (common::vterm::vt_beard77_v0(sqrt(rw2) * si::metres) / si::metres_per_second);
              break;

            case(vt_t::khvorostyanov_spherical):
              vt = common::vterm::vt_khvorostyanov(
                sqrt(rw2)           * si::metres, // TODO: consider caching rw?
                thrust::get<0>(tpl) * si::kelvins,
                thrust::get<2>(tpl) * si::kilograms / si::cubic_metres,
                thrust::get<3>(tpl) * si::pascals * si::seconds,
                true
              ) / si::metres_per_second;
              break;

            case(vt_t::khvorostyanov_nonspherical):
              vt = common::vterm::vt_khvorostyanov(
                sqrt(rw2)           * si::metres, // TODO: consider caching rw?
                thrust::get<0>(tpl) * si::kelvins,
                thrust::get<2>(tpl) * si::kilograms / si::cubic_metres,
                thrust::get<3>(tpl) * si::pascals * si::seconds,
                false
              ) / si::metres_per_second;
              break;

            default:
              vt = 0.; //sanity checks done in pimpl constructor
              break;
          }
          return vt;
        }   
      }; 

      template <typename real_t>
      struct common__vterm__vt__cached
      {
        vt_t vt_eq; //type of terminal velocity formula to use 

        //ctor
        common__vterm__vt__cached(const vt_t &vt_eq): vt_eq(vt_eq) {}

        BOOST_GPU_ENABLED 
        real_t operator()(
          const real_t &rw2, 
          const thrust::tuple<real_t, real_t, real_t, real_t> &tpl // (vt_0, p, rhod, eta)
        ) {   
#if !defined(__NVCC__)
          using std::sqrt;
#endif

          real_t vt;

          switch(vt_eq)
          {
            case(vt_t::beard77fast):
              vt = 
                common::vterm::vt_beard77_fact(
                  sqrt(rw2)           * si::metres,
                  thrust::get<1>(tpl) * si::pascals,
                  thrust::get<2>(tpl) * si::kilograms / si::cubic_metres,
                  thrust::get<3>(tpl) * si::pascals * si::seconds
                ) * thrust::get<0>(tpl); // cached vt_0
              break;

            default:
              vt = 0.; //sanity checks done in pimpl constructor
          }
          return vt;
        }   
      };

      template <typename real_t>
      struct common__vterm__ice
      {
        BOOST_GPU_ENABLED
        real_t operator()(
        const real_t &ice_a,
        const real_t &ice_c,
        const real_t &ice_rho,
         const thrust::tuple<real_t, real_t, real_t, real_t> &tpl // (T, p, rhod, eta)
        ) {
          return common::vterm::vt_beard76(
               ice_a * si::metres,
               thrust::get<0>(tpl) * si::kelvins,
               thrust::get<1>(tpl) * si::pascals,
               thrust::get<2>(tpl) * si::kilograms / si::cubic_metres,
               thrust::get<3>(tpl) * si::pascals * si::seconds
             ) / si::metres_per_second
            * real_t(common::moist_air::rho_i<real_t>() / common::moist_air::rho_w<real_t>());
      }
    };

      // if ice_switch=True, use this struct to select the correct formula
      template <typename real_t>
      struct select_vt
      {
        vt_t vt_eq;

        select_vt(const vt_t &vt_eq): vt_eq(vt_eq) {}

        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<
            real_t, real_t, real_t, real_t,  // rw2, ice_a, ice_c, rho_i
            real_t, real_t, real_t, real_t   // T, p, rhod, eta
          >& tpl) const
        {
          const real_t &rw2   = thrust::get<0>(tpl);
          const real_t &ice_a = thrust::get<1>(tpl);
          const real_t &ice_c = thrust::get<2>(tpl);
          const real_t &rho_i = thrust::get<3>(tpl);
          thrust::tuple<real_t, real_t, real_t, real_t> env(
            thrust::get<4>(tpl), thrust::get<5>(tpl), thrust::get<6>(tpl), thrust::get<7>(tpl)
          );

          if(rw2 == 0 && ice_a > 0 && ice_c > 0) {
            return common__vterm__ice<real_t>()(ice_a, ice_c, rho_i, env);
          }
          else if(vt_eq == vt_t::beard77fast) {
            return common__vterm__vt__cached<real_t>(vt_eq)(rw2, env);
          }
          else {
            return common__vterm__vt<real_t>(vt_eq)(rw2, env);
          }
        }
      };

   };


    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_vterm_invalid()
    {   
      namespace arg = thrust::placeholders;

      if(opts_init.terminal_velocity == vt_t::beard77fast) //use cached vt at sea level
      {
        thrust_device::vector<thrust_size_t> &vt0_bin(tmp_device_size_part);
        // get cached bin number
        thrust::transform_if(
          rw2.begin(), rw2.end(),
          vt.begin(),
          vt0_bin.begin(),
          detail::get_vt0_bin<real_t>(config.vt0_ln_r_min, config.vt0_ln_r_max, config.vt0_n_bin),
          arg::_1 == real_t(detail::invalid)
        );
        // calc the vt
        thrust::transform_if(
          rw2.begin(), rw2.end(),                                 // input - 1st arg
          thrust::make_zip_iterator(thrust::make_tuple(
            thrust::make_permutation_iterator(vt_0.begin(), vt0_bin.begin()),
            thrust::make_permutation_iterator(p.begin(),    ijk.begin()),
            thrust::make_permutation_iterator(rhod.begin(), ijk.begin()),
            thrust::make_permutation_iterator(eta.begin(),  ijk.begin())
          )),                                                     // input - 2nd arg   
          vt.begin(),                                             // condition argument
          vt.begin(),                                             // output
          detail::common__vterm__vt__cached<real_t>(opts_init.terminal_velocity),
          arg::_1 == real_t(detail::invalid)
        );
      }
      // non-cached vt
      else
        thrust::transform_if(
          rw2.begin(), rw2.end(),                                 // input - 1st arg
          thrust::make_zip_iterator(thrust::make_tuple(
            thrust::make_permutation_iterator(T.begin(),    ijk.begin()),
            thrust::make_permutation_iterator(p.begin(),    ijk.begin()),
            thrust::make_permutation_iterator(rhod.begin(), ijk.begin()),
            thrust::make_permutation_iterator(eta.begin(),  ijk.begin())
          )),                                                     // input - 2nd arg   
          vt.begin(),                                             // condition argument
          vt.begin(),                                             // output
          detail::common__vterm__vt<real_t>(opts_init.terminal_velocity),
          arg::_1 == real_t(detail::invalid)
        );
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_vterm_all()
    {
      if (opts_init.ice_switch)
      {
        auto zipped = thrust::make_zip_iterator(thrust::make_tuple(
            rw2.begin(),
            ice_a.begin(),
            ice_c.begin(),
            ice_rho.begin(),
            thrust::make_permutation_iterator(T.begin(),    ijk.begin()),
            thrust::make_permutation_iterator(p.begin(),    ijk.begin()),
            thrust::make_permutation_iterator(rhod.begin(), ijk.begin()),
            thrust::make_permutation_iterator(eta.begin(),  ijk.begin())
          ));
        auto zipped_end = zipped + rw2.size();

        thrust::transform(
          zipped,
          zipped_end,
          vt.begin(),
          detail::select_vt<real_t>(opts_init.terminal_velocity)
        );
      }
      else
      {
        if(opts_init.terminal_velocity == vt_t::beard77fast) //use cached vt at sea level
        {
          thrust_device::vector<thrust_size_t> &vt0_bin(tmp_device_size_part);
          // get cached bin number
          thrust::transform(
            rw2.begin(), rw2.end(),
            vt0_bin.begin(),
            detail::get_vt0_bin<real_t>(config.vt0_ln_r_min, config.vt0_ln_r_max, config.vt0_n_bin)
          );
          // calc the vt
          thrust::transform(
            rw2.begin(), rw2.end(),                                 // input - 1st arg
            thrust::make_zip_iterator(thrust::make_tuple(
              thrust::make_permutation_iterator(vt_0.begin(), vt0_bin.begin()),
              thrust::make_permutation_iterator(p.begin(),    ijk.begin()),
              thrust::make_permutation_iterator(rhod.begin(), ijk.begin()),
              thrust::make_permutation_iterator(eta.begin(),  ijk.begin())
            )),                                                     // input - 2nd arg
            vt.begin(),                                             // output
            detail::common__vterm__vt__cached<real_t>(opts_init.terminal_velocity)
          );
        }
        else
          thrust::transform(
            rw2.begin(), rw2.end(),                                 // input - 1st arg
            thrust::make_zip_iterator(thrust::make_tuple(
              thrust::make_permutation_iterator(T.begin(),    ijk.begin()),
              thrust::make_permutation_iterator(p.begin(),    ijk.begin()),
              thrust::make_permutation_iterator(rhod.begin(), ijk.begin()),
              thrust::make_permutation_iterator(eta.begin(),  ijk.begin())
            )),                                                     // input - 2nd arg
            vt.begin(),                                             // output
            detail::common__vterm__vt<real_t>(opts_init.terminal_velocity)
          );
      }
    }
  };  
};
