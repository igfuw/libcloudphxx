#pragma once

#include <libcloudph++/common/units.hpp>
#include <libcloudph++/common/macros.hpp>
#include <libcloudph++/common/earth.hpp>

// TODO: rename all such files to paper labels?

namespace libcloudphxx
{
  namespace common
  {
    namespace vterm
    {
      /// @brief dynamic viscosity of air (@copydetails Rogers_and_Yau_1989 third edition page 102) 
      /// temp in K and rhoa in kg/m3
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::dynamic_viscosity, real_t> visc(
	quantity<si::temperature, real_t> T
      ) {   
	return real_t(1.72 * 1e-5) * (real_t(393) / ( (T / si::kelvins) + real_t(120)) ) 
	  * real_t(pow(T/si::kelvins/273, real_t(3./2))) * si::kilograms / si::metres / si::seconds ;
      }

      // terminal fall velocity of spherical droplets 
      // TODO add another parametrisation for larger (nonspherical) drops
      // for derivation see @copydetails Khvorostyanov_and_Curry_2002 J. Atmos. Sci 
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::velocity, real_t> vt( 
	quantity<si::length, real_t> r, //radius
	quantity<si::temperature, real_t> T, //temperature
	quantity<si::mass_density, real_t> rhoa //density of air
      ) {
        using moist_air::rho_w;
        using earth::g;

	/// Best number (eq 2.7 in @copydetails Khvorostyanov_and_Curry_2002 J. Atmos. Sci) 
	/// with maximum projected cross-sectional area parametrised as for spherical droplets (A=pi/4*D^2)
	quantity<si::dimensionless, real_t> X = real_t(32./3) * (rho_w<real_t>() - rhoa)/rhoa * g<real_t>() 
	  * r * r * r / visc(T) / visc(T) * rhoa * rhoa; //TODO use pow<>()  

	/// terminal velocity parametrisation coeffs 
	/// eqs 2.12, 2.13 in @copydetails Khvorostyanov_and_Curry_2002 J. Atmos. Sci 
	quantity<si::dimensionless, real_t> b = real_t(.0902/2) * sqrt(X)
	  * pow(sqrt(real_t(1)+real_t(.0902)*sqrt(X))-real_t(1), -1)
	  * pow(sqrt(real_t(1)+real_t(.0902)*sqrt(X)), -1) ;
	quantity<si::dimensionless, real_t> a = real_t(9.06 * 9.06 / 4)
	  * pow(sqrt(real_t(1)+real_t(.0902)*sqrt(X))-real_t(1), 2) / pow(X,b) ;

	/// eq 3.1 in @copydetails Khvorostyanov_and_Curry_2002 J. Atmos. Sci 
	quantity<si::dimensionless, real_t> Av = a
	  * pow(visc(T) / rhoa * real_t(1e4) * si::seconds / si::square_metres, real_t(1)-real_t(2)*b)
	  * pow(real_t(4./3) * rho_w<real_t>() / rhoa * g<real_t>() *real_t(1e2)* si::seconds * si::seconds / si::metres, b) ;
	quantity<si::dimensionless, real_t> Bv = real_t(3)*b - real_t(1) ;

	return (Av * real_t(pow(real_t(2*1e2) * r/si::metres, Bv)))/real_t(1e2) * si::metres_per_second  ;
      }
    };
  };
};
