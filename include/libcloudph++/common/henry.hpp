#pragma once

#include <libcloudph++/common/units.hpp>
#include <libcloudph++/common/macros.hpp> 
#include <libcloudph++/common/earth.hpp> 
#include <libcloudph++/common/molar_mass.hpp>
#include <libcloudph++/common/moist_air.hpp>

namespace libcloudphxx
{
  namespace common
  {
    typedef divide_typeof_helper<
      divide_typeof_helper<
        si::amount,  
        si::volume
      >::type,
      si::pressure
    >::type amount_over_volume_over_pressure;

   typedef divide_typeof_helper<
      si::dimensionless,  
      si::time
    >::type one_over_time;

    namespace henry
    {
      // Henry constant [1e3 (per litre); p_stp (per standard atmosphere)]
      libcloudphxx_const_derived(amount_over_volume_over_pressure, H_SO2,  real_t(1.23        * 1e3) / earth::p_stp<real_t>() * si::moles / si::cubic_metres)
      libcloudphxx_const_derived(amount_over_volume_over_pressure, H_H2O2, real_t(7.45e4      * 1e3) / earth::p_stp<real_t>() * si::moles / si::cubic_metres)
      libcloudphxx_const_derived(amount_over_volume_over_pressure, H_O3,   real_t(1.13 * 1e-2 * 1e3) / earth::p_stp<real_t>() * si::moles / si::cubic_metres)
      libcloudphxx_const_derived(amount_over_volume_over_pressure, H_NH3,  real_t(62          * 1e3) / earth::p_stp<real_t>() * si::moles / si::cubic_metres)
      libcloudphxx_const_derived(amount_over_volume_over_pressure, H_HNO3, real_t(2.1 * 1e5   * 1e3) / earth::p_stp<real_t>() * si::moles / si::cubic_metres)
      libcloudphxx_const_derived(amount_over_volume_over_pressure, H_CO2,  real_t(3.4 * 1e-2  * 1e3) / earth::p_stp<real_t>() * si::moles / si::cubic_metres)

      //modification to Henry due to temperature - Table 4 in Kreidenweiss et al 2003
      libcloudphxx_const_derived(si::temperature, dHR_SO2,  real_t(3150.) * si::kelvins)
      libcloudphxx_const_derived(si::temperature, dHR_O3,   real_t(2540.) * si::kelvins)
      libcloudphxx_const_derived(si::temperature, dHR_H2O2, real_t(7300.) * si::kelvins)
      libcloudphxx_const_derived(si::temperature, dHR_NH3,  real_t(4100.) * si::kelvins)
      libcloudphxx_const_derived(si::temperature, dHR_HNO3, real_t(8700.) * si::kelvins)
      libcloudphxx_const_derived(si::temperature, dHR_CO2,  real_t(2440.) * si::kelvins)

      //gas phase diffusivity - 
      //p1: W. J. Massman, Atmospheric Environment Vol. 32, No. 6, pp. 1111—1127, 1998
      //    A review of the molecular diffusivities of H O, CO , CH , CO, O , SO , NH , N O, NO, and NO in air, O and N near stp
      //    Table 2
      //p2: M. J. Tang et al, Atmos. Chem. Phys., 14, 9233–9247, 2014
      //    Compilation and evaluation of gas phase diffusion coefficients of reactive trace gases 
      //    in the atmosphere: volume 1. Inorganic compounds
      //    Table 2
      libcloudphxx_const_derived(diffusivity, D_SO2,  real_t(.1089 * 1e-4) * si::metres * si::metres / si::seconds) //p1
      libcloudphxx_const_derived(diffusivity, D_O3,   real_t(.1444 * 1e-4) * si::metres * si::metres / si::seconds) //p1
      libcloudphxx_const_derived(diffusivity, D_H2O2, real_t(.8700 * 1e-4) * si::metres * si::metres / si::seconds) //p2
      libcloudphxx_const_derived(diffusivity, D_CO2,  real_t(.1381 * 1e-4) * si::metres * si::metres / si::seconds) //p1
      libcloudphxx_const_derived(diffusivity, D_HNO3, real_t(.6525 * 1e-4) * si::metres * si::metres / si::seconds) //p2
      libcloudphxx_const_derived(diffusivity, D_NH3,  real_t(.1978 * 1e-4) * si::metres * si::metres / si::seconds) //p1

      //accomodation coefficient - Table 4 in Kreidenweiss et al 2003
      libcloudphxx_const_derived(si::dimensionless, ac_SO2,  real_t(.035)  )
      libcloudphxx_const_derived(si::dimensionless, ac_O3,   real_t(.00053))
      libcloudphxx_const_derived(si::dimensionless, ac_H2O2, real_t(.018)  )
      libcloudphxx_const_derived(si::dimensionless, ac_CO2,  real_t(.05)   )
      libcloudphxx_const_derived(si::dimensionless, ac_HNO3, real_t(.05)   )
      libcloudphxx_const_derived(si::dimensionless, ac_NH3,  real_t(.05)   )

      // mean molecular velocity
      // mean of the Maxwell distribution at a given temperature
      // eq. 8.2 from Seinfeld and Pandis
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::velocity, real_t> molec_vel(
        const quantity<si::temperature, real_t> &T,
        const quantity<mass_over_amount, real_t> &M
      ) {
        return (
          real_t(  //bug in boost #6957
            sqrt(
              real_t(8.) / 
#if !defined(__NVCC__)
	      pi<real_t>()
#else
	      CUDART_PI
#endif
              * (moist_air::kaBoNA<real_t>() * T  / M) / si::square_metres * si::seconds * si::seconds
            )
          ) * si::metres / si::seconds
        );
      }

      //mass transfer coefficient
      //Peter Warneck - Chemistry of the Natural Atmosphere (chapter 8.4.2 eq. 8.23)
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<one_over_time, real_t> mass_trans(
        const quantity<si::area, real_t> &rw2,
        const quantity<diffusivity, real_t> &D,
        const quantity<si::dimensionless, real_t> &acc_coeff,
        const quantity<si::temperature, real_t> &T,
        const quantity<mass_over_amount, real_t> &M
      ) {
        return(
          real_t(1.) / (
            rw2 / real_t(3.) / D 
            + 
            real_t(4./3.) / acc_coeff * (real_t(sqrt(rw2/si::square_meter)) * si::metres) / molec_vel(T, M))
        );                               //bug in boost #6957
      }

      //Henry constant(temperature)
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<amount_over_volume_over_pressure, real_t> H_temp(
        const quantity<si::temperature, real_t> &T,
        const quantity<amount_over_volume_over_pressure, real_t> &H,
        const quantity<si::temperature, real_t> &dHR
      ) {
        return (H * exp(dHR * (real_t(1.) / T - (real_t(1./298) / si::kelvins)))); 
      }
    };
  };
};
