#include <libcloudph++/common/moist_air.hpp>

#if defined(__NVCC__)
#  include <math_constants.h>
#endif

namespace libcloudphxx
{
  namespace lgrngn
  {
    using detail::tpl_calc_wrap;

    template <typename real_t, typename n_t>
    struct kernel_base
    {
      // pointer to kernel parameters device vector
      thrust_device::pointer<real_t> k_params;

      // number of user-defined parameters
      n_t n_user_params; 

      // largest radius for which efficiency is defined, 0 - n/a
      real_t r_max;
   
      //ctor
      BOOST_GPU_ENABLED
      kernel_base(thrust_device::pointer<real_t> k_params, n_t n_user_params = 0, real_t r_max = 0.) : 
        k_params(k_params), n_user_params(n_user_params), r_max(r_max) {}

      // thrust requires that a default ctor exists
      BOOST_GPU_ENABLED
      kernel_base() = default;

      BOOST_GPU_ENABLED
      virtual real_t calc(const tpl_calc_wrap<real_t,n_t> &) const {return 0;}
    };


    //Golovin kernel
    template <typename real_t, typename n_t>
    struct kernel_golovin : kernel_base<real_t, n_t>
    {
      //ctor
      BOOST_GPU_ENABLED
      kernel_golovin(thrust_device::pointer<real_t> k_params) : kernel_base<real_t, n_t>(k_params, 1) {}

      // thrust requires that a default ctor exists
      BOOST_GPU_ENABLED
      kernel_golovin() = default;

      BOOST_GPU_ENABLED
      virtual real_t calc(const tpl_calc_wrap<real_t,n_t> &tpl_wrap) const
      {
        enum { n_a_ix, n_b_ix, rw2_a_ix, rw2_b_ix, vt_a_ix, vt_b_ix, rd3_a_ix, rd3_b_ix };
#if !defined(__NVCC__)
        using std::max;
        using std::sqrt;
#endif

        const real_t rw2_a = thrust::get<rw2_a_ix>(tpl_wrap.get_rw());
        const real_t rw2_b = thrust::get<rw2_b_ix>(tpl_wrap.get_rw());

        real_t res =
#if !defined(__NVCC__)
        pi<real_t>()
#else
        CUDART_PI
#endif
        * 4. / 3.
        * kernel_base<real_t, n_t>::k_params[0]
        * max(
            thrust::get<n_a_ix>(tpl_wrap.get_rw()),
            thrust::get<n_b_ix>(tpl_wrap.get_rw())
          )
        * (
            rw2_a * sqrt(rw2_a) +
            rw2_b * sqrt(rw2_b) 
          );
        return res;
      }
    };


    //geometric kernel
    template <typename real_t, typename n_t>
    struct kernel_geometric : kernel_base<real_t, n_t>
    {
      //ctor (default one)
      BOOST_GPU_ENABLED
      kernel_geometric(thrust_device::pointer<real_t> k_params = thrust_device::pointer<real_t>(), n_t n_user_params = 0, real_t r_max = 0.) : 
        kernel_base<real_t, n_t>(k_params, n_user_params, r_max) {}

      //bilinear interpolation of collision efficiencies, required by dervied classes
      BOOST_GPU_ENABLED
      real_t interpolated_efficiency(real_t, real_t) const;

      BOOST_GPU_ENABLED
      virtual real_t calc(const tpl_calc_wrap<real_t,n_t> &tpl_wrap) const
      {
        enum { n_a_ix, n_b_ix, rw2_a_ix, rw2_b_ix, vt_a_ix, vt_b_ix, rd3_a_ix, rd3_b_ix };
#if !defined(__NVCC__)
        using std::abs;
        using std::max;
        using std::sqrt;
#endif
        return 
#if !defined(__NVCC__)
        pi<real_t>()
#else
        CUDART_PI
#endif
        * max(
            thrust::get<n_a_ix>(tpl_wrap.get_rw()),
            thrust::get<n_b_ix>(tpl_wrap.get_rw())
          )
        * abs(
            thrust::get<vt_a_ix>(tpl_wrap.get_rw()) -
            thrust::get<vt_b_ix>(tpl_wrap.get_rw())
          )
        * (thrust::get<rw2_a_ix>(tpl_wrap.get_rw()) +
           thrust::get<rw2_b_ix>(tpl_wrap.get_rw()) +
           2.*sqrt(thrust::get<rw2_a_ix>(tpl_wrap.get_rw())*thrust::get<rw2_b_ix>(tpl_wrap.get_rw()))
          );
    //    return res;
      }
    };

    //geometric kernel with a multiplier
    template <typename real_t, typename n_t>
    struct kernel_geometric_with_multiplier : kernel_geometric<real_t, n_t>
    {
      //ctor
      BOOST_GPU_ENABLED
      kernel_geometric_with_multiplier(thrust_device::pointer<real_t> k_params) : kernel_geometric<real_t, n_t>(k_params, 1) {}

      // thrust requires that a default ctor exists
      BOOST_GPU_ENABLED
      kernel_geometric_with_multiplier() = default;

      BOOST_GPU_ENABLED
      virtual real_t calc(const tpl_calc_wrap<real_t,n_t> &tpl_wrap) const
      {
        return kernel_geometric<real_t, n_t>::calc(tpl_wrap) * kernel_base<real_t, n_t>::k_params[0];
      }
    };

    //Long kernel
    template <typename real_t, typename n_t>
    struct kernel_long : kernel_geometric<real_t, n_t>
    {
      enum { n_a_ix, n_b_ix, rw2_a_ix, rw2_b_ix, vt_a_ix, vt_b_ix, rd3_a_ix, rd3_b_ix };
      //ctor
      BOOST_GPU_ENABLED
      kernel_long() : kernel_geometric<real_t, n_t>() {}

      BOOST_GPU_ENABLED
      virtual real_t calc(const tpl_calc_wrap<real_t,n_t> &tpl_wrap) const
      {
#if !defined(__NVCC__)
        using std::abs;
        using std::max;
        using std::min;
        using std::sqrt;
#endif
        real_t res = kernel_geometric<real_t, n_t>::calc(tpl_wrap);

        real_t r_L = max(sqrt(thrust::get<rw2_a_ix>(tpl_wrap.get_rw())), sqrt(thrust::get<rw2_b_ix>(tpl_wrap.get_rw())));
        if(r_L < 50.e-6)
        {
          real_t r_s = min(sqrt(thrust::get<rw2_a_ix>(tpl_wrap.get_rw())), sqrt(thrust::get<rw2_b_ix>(tpl_wrap.get_rw())));
          if(r_s <= 3e-6)
            res = 0.;
          else
            res *= 4.5e8 * r_L * r_L * (1. - 3e-6/r_s);
        }
        
        return  res;
      }
    };

    template <typename real_t, typename n_t>
    struct kernel_geometric_with_efficiencies : kernel_geometric<real_t, n_t>
    {
      //ctor
      BOOST_GPU_ENABLED
      kernel_geometric_with_efficiencies(thrust_device::pointer<real_t> k_params, real_t r_max) : kernel_geometric<real_t, n_t>(k_params, 0, r_max) {}

      // thrust requires that a default ctor exists
      BOOST_GPU_ENABLED
      kernel_geometric_with_efficiencies() = default;

      BOOST_GPU_ENABLED
      virtual real_t calc(const tpl_calc_wrap<real_t,n_t> &tpl_wrap) const
      {
        enum { n_a_ix, n_b_ix, rw2_a_ix, rw2_b_ix, vt_a_ix, vt_b_ix, rd3_a_ix, rd3_b_ix };

#if !defined(__NVCC__)
        using std::sqrt;
#endif

        return  kernel_geometric<real_t, n_t>::interpolated_efficiency(
                  sqrt( thrust::get<rw2_a_ix>(tpl_wrap.get_rw())),
                  sqrt( thrust::get<rw2_b_ix>(tpl_wrap.get_rw()))
                ) * kernel_geometric<real_t, n_t>::calc(tpl_wrap);
      }
    };

    // turbulent kernel from the 2015 JAS Onishi paper
    // two user params defined at initialization: 
    // turbulence dissipataion rate and Taylor microscale Reynolds number
    // TODO: get these values from flow characteristic (simulation-time during hskpng)
    //       cf. Benmoshe et al, JGR 2012
    template <typename real_t, typename n_t>
    struct kernel_onishi : kernel_geometric<real_t, n_t>
    {
      detail::wang_collision_enhancement_t<real_t> wang_collision_enhancement;

      //ctor
      BOOST_GPU_ENABLED
      kernel_onishi(thrust_device::pointer<real_t> k_params, real_t r_max) : kernel_geometric<real_t, n_t>(k_params, 1, r_max) {}

      // thrust requires that a default ctor exists
      BOOST_GPU_ENABLED
      kernel_onishi() = default;

      BOOST_GPU_ENABLED
      virtual real_t calc(const tpl_calc_wrap<real_t,n_t> &tpl_wrap) const
      {
        enum { n_a_ix, n_b_ix, rw2_a_ix, rw2_b_ix, vt_a_ix, vt_b_ix, rd3_a_ix, rd3_b_ix };
        enum { rhod_ix, eta_ix, diss_rate_ix };

#if !defined(__NVCC__)
        using std::sqrt;
#endif
        real_t rwa = sqrt( thrust::get<rw2_a_ix>(tpl_wrap.get_rw()));
        real_t rwb = sqrt( thrust::get<rw2_b_ix>(tpl_wrap.get_rw()));
        real_t onishi_nograv = detail::kernel_onishi_nograv<real_t>           // value of the turbulent onishi kernel that neglects gravitational settling
        (
          rwa, rwb, kernel_base<real_t, n_t>::k_params[0], thrust::get<diss_rate_ix>(tpl_wrap.get_ro_calc()),                  // k_params[0] - Re_lambda
          thrust::get<eta_ix>(tpl_wrap.get_ro_calc()) / thrust::get<rhod_ix>(tpl_wrap.get_ro_calc()),                          // kinetic viscosity 
          common::moist_air::rho_w<real_t>() / si::kilograms * si::cubic_metres / thrust::get<rhod_ix>(tpl_wrap.get_ro_calc()) // ratio of water to air density
        );

        real_t geometric = kernel_geometric<real_t, n_t>::calc(tpl_wrap);
        real_t res = 
          kernel_geometric<real_t, n_t>::interpolated_efficiency(rwa, rwb) *             // stagnant air collision efficiency
          wang_collision_enhancement(rwa, rwb, kernel_base<real_t, n_t>::k_params[0]) *  // Wang turbulent collision efficiency enhancement, k_params[0] - epsilon
          sqrt(
            geometric * geometric +  // geometric kernel 
            onishi_nograv * onishi_nograv
          );

        return res;
      }
    };
  };
};


                            
