#pragma once

#include <libcloudph++/common/units.hpp>
#include <libcloudph++/common/macros.hpp>
#include <libcloudph++/common/earth.hpp>
#include <libcloudph++/common/const_cp.hpp>

// TODO: rename all such files to paper labels?

namespace libcloudphxx
{
  namespace common
  {
    namespace vterm
    {
// TODO: move visc to a separate file
      /// @brief dynamic viscosity of air (@copydetails Rogers_and_Yau_1989 third edition page 102) 
      /// temp in K and rhoa in kg/m3
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::dynamic_viscosity, real_t> visc(
	quantity<si::temperature, real_t> T
      ) 
      {  
        using const_cp::T_tri;
 
	return real_t(1.72 * 1e-5) * (real_t(393) / ( (T / si::kelvins) + real_t(120)) ) 
	  * real_t(pow(T/T_tri<real_t>(), real_t(3./2))) * si::kilograms / si::metres / si::seconds;
      }

      // terminal fall velocity of spherical droplets 
      // TODO add another parametrisation for larger (nonspherical) drops
      // for derivation see @copydetails Khvorostyanov_and_Curry_2002 J. Atmos. Sci 
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::velocity, real_t> vt( 
	quantity<si::length, real_t> r, //radius
	quantity<si::temperature, real_t> T, //temperature
	quantity<si::mass_density, real_t> rhoa, //density of air
        quantity<si::dynamic_viscosity, real_t> eta
      ) {
        using moist_air::rho_w;
        using earth::g;

	/// Best number (eq 2.7 in @copydetails Khvorostyanov_and_Curry_2002 J. Atmos. Sci) 
	/// with maximum projected cross-sectional area parametrised as for spherical droplets (A=pi/4*D^2)
	quantity<si::dimensionless, real_t> X = real_t(32./3) * (rho_w<real_t>() - rhoa)/rhoa * g<real_t>() 
	  * r * r * r / eta / eta * rhoa * rhoa; //TODO use pow<>()  

	/// terminal velocity parametrisation coeffs 
	/// eqs 2.12, 2.13 in @copydetails Khvorostyanov_and_Curry_2002 J. Atmos. Sci 
	quantity<si::dimensionless, real_t> b = real_t(.0902/2) * sqrt(X)
	  * pow(sqrt(real_t(1)+real_t(.0902)*sqrt(X))-real_t(1), -1)
	  * pow(sqrt(real_t(1)+real_t(.0902)*sqrt(X)), -1) ;
	quantity<si::dimensionless, real_t> a = real_t(9.06 * 9.06 / 4)
	  * pow(sqrt(real_t(1)+real_t(.0902)*sqrt(X))-real_t(1), 2) / pow(X,b) ;

	/// eq 3.1 in @copydetails Khvorostyanov_and_Curry_2002 J. Atmos. Sci 
	quantity<si::dimensionless, real_t> Av = a
	  * pow(eta / rhoa * real_t(1e4) * si::seconds / si::square_metres, real_t(1)-real_t(2)*b)
	  * pow(real_t(4./3) * rho_w<real_t>() / rhoa * g<real_t>() *real_t(1e2)* si::seconds * si::seconds / si::metres, b) ;
	quantity<si::dimensionless, real_t> Bv = real_t(3)*b - real_t(1) ;

	return (Av * real_t(pow(real_t(2*1e2) * r/si::metres, Bv)))/real_t(1e2) * si::metres_per_second  ;
      }

      // terminal fall velocity according to Beard (1976)
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::velocity, real_t> vt_beard( 
	quantity<si::length, real_t> r, //radius
	quantity<si::temperature, real_t> T, //temperature
	quantity<si::pressure, real_t> p, //pressure
	quantity<si::mass_density, real_t> rhoa, //density of air
	quantity<si::mass_density, real_t> rhow, //density of water
        quantity<si::dynamic_viscosity, real_t> eta
      ) 
      {
        if(r <= quantity<si::length, real_t>(9.5e-6 * si::meter)) //TODO: < 0.5um
        {
          using earth::p_stp;
          using earth::g;
          quantity<si::length, real_t> l = 6.62e-8 * si::meter * (eta / 1.818e-5 / si::pascal / si::seconds) * (p_stp<real_t>() / p) * sqrt(T / 293.15 / si::kelvin);
          quantity<si::dimensionless, real_t> C_ac = 1. + 1.255 * l / r;
          return ( (rhow<real_t>()-rhoa<real_t>()) * g<real_t>() / ( 18. * eta) * C_ac * r);
        } 
        else if(r <= quantity<si::length, real_t>(5.035e-4 * si::meter))
        {
          const real_t b[7] = { -0.318657e1, 0.992696, -0.153193e-2, -0.987059e-3, -0.578878e-3, 0.855176e-4,-0.327815e-5};

        }
          const real_t cdata c /-0.500015e1,0.523778e1,-0.204914e1,0.475294,-0.542819e-1,
      &         0.238449e-2/
         
      }
    
    };
  };
};
