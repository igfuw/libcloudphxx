#pragma once

#include "units.hpp"
#include "macros.hpp"

namespace libcloudphxx
{
  namespace common
  {
    namespace ventil
    {
      // Reynolds number for a particles falling with terminal velocity
      // see e.g. section 4 in Smolik et al 2001, Aerosol Sci.
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::dimensionless, real_t> Re(
	const quantity<si::velocity, real_t> v_term,       // particle terminal velocity
        const quantity<si::length, real_t> r_w,            // particle wet radius
        const quantity<si::mass_density, real_t> rho,      // air density
        const quantity<si::dynamic_viscosity, real_t>  eta // air viscosity 
      )
      {
	return v_term * (real_t(2) * r_w) * rho / eta;
      }

      // Nusselt number
      // see eq. 1 in Smolik et al 2001, Aerosol Sci.
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::dimensionless, real_t> Nu(
        const quantity<si::dimensionless, real_t> Pr, // Prandtl number
        const quantity<si::dimensionless, real_t> Re  // Reynolds number
      ) 
      {
#if !defined(__NVCC__)
	//using std::pow;
	using std::max;
#endif

        return real_t(1) + cbrt(real_t(1) + Re * Pr) * max(real_t(1), pow(real_t(Re), real_t(.077)));
        //                                  ^^^^^^^ 
        //                      Peclet number /
      }

      // Sherwood number
      // see eq. 2 in Smolik et al 2001, Aerosol Sci.
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::dimensionless, real_t> Sh(
        const quantity<si::dimensionless, real_t> Sc, // Schmidt number
        const quantity<si::dimensionless, real_t> Re  // Reynolds number
      ) 
      {
        return Nu(Sc, Re);
      }

      // Schmidt number
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::dimensionless, real_t> Sc(
        const quantity<si::dynamic_viscosity, real_t> &eta,
        const quantity<si::mass_density, real_t> &rho,
        const quantity<diffusivity, real_t> &D
      )
      {
        return eta / rho / D;
      }

      // Prandtl number
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::dimensionless, real_t> Pr(
        const quantity<si::dynamic_viscosity, real_t> &eta,
        const quantity<energy_over_temperature_mass, real_t> &c_p,
        const quantity<thermal_conductivity, real_t> &K
      )
      {
        return c_p * eta / K;
      }
    };
  };
};
