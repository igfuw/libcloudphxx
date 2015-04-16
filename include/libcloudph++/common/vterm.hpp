#pragma once

#include <libcloudph++/common/units.hpp>
#include <libcloudph++/common/macros.hpp>
#include <libcloudph++/common/earth.hpp>
#include <libcloudph++/common/const_cp.hpp>
#include <libcloudph++/common/kelvin_term.hpp>

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
      quantity<si::velocity, real_t> vt_khvorostyanov_spherical( 
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

      // terminal fall velocity according to Beard (1976), Table 1
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::velocity, real_t> vt_beard( 
	quantity<si::length, real_t> r, //radius
	quantity<si::temperature, real_t> T, //temperature
	quantity<si::pressure, real_t> p, //pressure
	quantity<si::mass_density, real_t> rhoa, //density of air
        quantity<si::dynamic_viscosity, real_t> eta
      ) 
      {
        using earth::p_stp;
        using earth::g;
        using moist_air::rho_w;

        if(r <= quantity<si::length, real_t>(real_t(9.5e-6) * si::meters)) //TODO: < 0.5um
        {
          quantity<si::dimensionless, real_t> l = ( real_t(6.62e-8)  * (eta / si::pascals / si::seconds/ real_t(1.818e-5) )  * (p_stp<real_t>() / p)  *  pow(T / si::kelvins / real_t(293.15), real_t(1./2.)) );
          quantity<si::dimensionless, real_t> C_ac = real_t(1.) + real_t(1.255) * l * si::meters / r;
          return ( (rho_w<real_t>()-rhoa) * g<real_t>() / ( real_t(4.5) * eta) * C_ac * r *r);
        } 

        else if(r <= quantity<si::length, real_t>(real_t(5.035e-4) * si::meters))
        {
          const real_t b[7] = { -0.318657e1, 0.992696, -0.153193e-2, -0.987059e-3, -0.578878e-3, 0.855176e-4,-0.327815e-5};
          quantity<si::dimensionless, real_t> l = ( real_t(6.62e-8)  * (eta / si::pascals / si::seconds/ real_t(1.818e-5) )  * (p_stp<real_t>() / p)  *  pow(T / si::kelvins / real_t(293.15), real_t(1./2.)) );
          quantity<si::dimensionless, real_t> C_ac = real_t(1.) + real_t(1.255) * l * si::meters / r;
          quantity<si::dimensionless, real_t> log_N_Da = log( real_t(32./3.) * r * r * r * rhoa * (rho_w<real_t>() - rhoa) * g<real_t>() / eta / eta );
          quantity<si::dimensionless, real_t> Y = 0.;
          for(int i=0; i<7; ++i)
            Y = Y + b[i] * pow(log_N_Da, real_t(i));
          quantity<si::dimensionless, real_t> N_Re = C_ac * exp(Y);
          return (eta * N_Re / rhoa / real_t(2.) / r);
        }
        else //TODO: > 7mm
        {
          using kelvin::sg_surf; 
          const real_t b[6] = { -0.500015e1, 0.523778e1, -0.204914e1, 0.475294, -0.542819e-1, 0.238449e-2};
          quantity<si::dimensionless, real_t> Bo = real_t(16./3.) * r * r * (rho_w<real_t>() - rhoa) * g<real_t>() / sg_surf<real_t>(T);
          quantity<si::dimensionless, real_t> N_p = sg_surf<real_t>(T) * sg_surf<real_t>(T) * sg_surf<real_t>(T)  * rhoa * rhoa / eta / eta / eta / eta / g<real_t>() / (rho_w<real_t>() - rhoa);
          quantity<si::dimensionless, real_t> X = log (Bo * pow(N_p, real_t(1./6.)));
          quantity<si::dimensionless, real_t> Y = 0.;
          for(int i=0; i<6; ++i)
            Y = Y + b[i] * pow(X, real_t(i));
          quantity<si::dimensionless, real_t> N_Re = pow(N_p, real_t(1./6.)) * exp(Y);
          return (eta * N_Re / rhoa / real_t(2.) / r);
        }         
      }    
    };
  };
};
