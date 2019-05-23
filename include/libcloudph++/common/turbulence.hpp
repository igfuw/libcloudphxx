#pragma once

namespace libcloudphxx
{
  namespace common
  {
    // Grabowski and Abade 2017
    namespace turbulence
    {
      libcloudphxx_const_derived(si::dimensionless, C_E, real_t(0.845));
      libcloudphxx_const_derived(si::dimensionless, C_tau, real_t(1.5));
      libcloudphxx_const_derived(si::dimensionless, cube_root_of_two_pi, pow(real_t(2) *
#if !defined(__NVCC__)
        pi<real_t>()
#else
        CUDART_PI
#endif
        , real_t(1./3.))
      );

#if !defined(__NVCC__)
      using std::pow;
      using std::sqrt;
#endif

      typedef divide_typeof_helper<
        si::dimensionless,  
        si::length
      >::type one_over_length;
      libcloudphxx_const(one_over_length, a_1, real_t(3e-4), real_t(1) / si::meters);

      typedef divide_typeof_helper<
        si::area,  
        si::time
       >::type area_over_time;
      libcloudphxx_const(area_over_time, a_2, real_t(2.8e-4), si::square_meters / si::seconds);

      typedef divide_typeof_helper<
        si::area,  
        multiply_typeof_helper<si::time, si::time>::type
      >::type area_over_time_squared;

      typedef divide_typeof_helper<
        area_over_time_squared, 
        si::time
      >::type area_over_time_cubed;

      typedef divide_typeof_helper<
        si::dimensionless, 
        si::time
      >::type one_over_time;

      typedef divide_typeof_helper<
        si::dimensionless, 
        si::area
      >::type one_over_area;

      template <typename real_t>
      quantity<si::length, real_t> length_scale(
        quantity<si::length, real_t> dx
      )
      { return quantity<si::length, real_t>(dx);}

      template <typename real_t>
      quantity<si::length, real_t> length_scale(
        quantity<si::length, real_t> dx,
        quantity<si::length, real_t> dz
      )
      //{ return quantity<si::length, real_t>(sqrt(dx*dz));}
      { return quantity<si::length, real_t>(dz);}

      template <typename real_t>
      quantity<si::length, real_t> length_scale(
        quantity<si::length, real_t> dx,
        quantity<si::length, real_t> dy,
        quantity<si::length, real_t> dz
      )
      //{ return quantity<si::length, real_t>(pow((dx/si::metres)*(dy/si::metres)*(dz/si::metres), real_t(1./3.)) * si::metres);}
      { return quantity<si::length, real_t>(dz);}

      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<area_over_time_squared, real_t> tke(
        const quantity<area_over_time_cubed, real_t> &diss_rate, // dissipation rate (epsilon)
        const quantity<si::length, real_t> &L // characteristic length-scale
      )
      {
        return quantity<area_over_time_squared, real_t>(pow((L * diss_rate) / si::cubic_metres * si::seconds * si::seconds * si::seconds / C_E<real_t>(), real_t(2./3.)) * si::metres * si::metres / si::seconds / si::seconds);
      };

      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::time, real_t> tau(
        const quantity<area_over_time_squared, real_t> &tke,
        const quantity<si::length, real_t> &L // characteristic length-scale
      )
      {
        return quantity<si::time, real_t>(L / cube_root_of_two_pi<real_t>()  * sqrt(C_tau<real_t>() / tke));
      };

      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::velocity, real_t> update_turb_vel(
        const quantity<si::velocity, real_t> &wp,
        const quantity<si::time, real_t> &tau,
        const quantity<si::time, real_t> &dt,
        const quantity<area_over_time_squared, real_t> &tke,
        const quantity<si::dimensionless, real_t> &r_normal
      )
      {
        quantity<si::dimensionless, real_t> exp_m_dt_ov_tau(exp(-dt / tau));
        return quantity<si::velocity, real_t>(wp * exp_m_dt_ov_tau + sqrt((real_t(1) - exp_m_dt_ov_tau * exp_m_dt_ov_tau) * real_t(2./3.) * tke ) * r_normal);
      };

      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::time, real_t> tau_relax(
        const quantity<one_over_area, real_t> &wet_mom_1_over_vol
      )
      {
        return quantity<si::time, real_t>(real_t(1) / ( a_2<real_t>() * wet_mom_1_over_vol ) ); 
      };

      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<one_over_time, real_t> dot_turb_ss(
        const quantity<si::dimensionless, real_t> &ssp,
        const quantity<si::velocity, real_t> &wp,
        const quantity<si::time, real_t> &tau_rlx
      )
      {
        return quantity<one_over_time, real_t>(a_1<real_t>() *  wp - ssp / tau_rlx);
      };
    };
  };
};
