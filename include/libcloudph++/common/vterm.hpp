#pragma once

#include "units.hpp"
#include "macros.hpp"
#include "earth.hpp"
#include "const_cp.hpp"
#include "kelvin_term.hpp"

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
        const real_t T_over_T_tri = T/T_tri<real_t>();
 
        return real_t(1.72 * 1e-5) * (real_t(393) / ( (T / si::kelvins) + real_t(120)) ) 
          * real_t(T_over_T_tri * sqrt(T_over_T_tri)) * si::kilograms / si::metres / si::seconds;
      }

      // terminal fall velocity of spherical droplets 
      // TODO add another parametrisation for larger (nonspherical) drops
      // for derivation see @copydetails Khvorostyanov_and_Curry_2002 J. Atmos. Sci 
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::velocity, real_t> vt_khvorostyanov( 
        quantity<si::length, real_t> r, //radius
        quantity<si::temperature, real_t> T, //temperature
        quantity<si::mass_density, real_t> rhoa, //density of air
        quantity<si::dynamic_viscosity, real_t> eta,
        bool spherical // flag if we include nonsphericity (0) or not (1)
      ) {
#if !defined(__NVCC__)
        using std::pow;
        using std::sqrt;
#endif
        using moist_air::rho_w;
        using earth::g;

        // explicitly convert to doubles
        quantity<si::length, double> r_dbl(r);//.value * si::metres); 
        quantity<si::mass_density, double> rhoa_dbl(rhoa);//.value * si::kilograms / si::cubic_metres); 
        quantity<si::dynamic_viscosity, double> eta_dbl(eta);//.value * si::pascal * si::second); 

        /// Best number (eq 2.7 in @copydetails Khvorostyanov_and_Curry_2002 J. Atmos. Sci) 
        /// with maximum projected cross-sectional area parametrised as for spherical droplets (A=pi/4*D^2)
        quantity<si::dimensionless, double> X = double(32./3) * (rho_w<double>() - rhoa_dbl)/rhoa_dbl * g<double>() 
          * r_dbl * r_dbl * r_dbl / eta_dbl / eta_dbl * rhoa_dbl * rhoa_dbl; //TODO use pow<>()  

        /// terminal velocity parametrisation coeffs 
        /// eqs 2.12, 2.13 in @copydetails Khvorostyanov_and_Curry_2002 J. Atmos. Sci 
        quantity<si::dimensionless, double> b = double(.0902/2) * sqrt(X) / 
          ( (sqrt(double(1)+double(.0902)*sqrt(X))-double(1))
          * (sqrt(double(1)+double(.0902)*sqrt(X)))) ;
        const double pow_hlpr = sqrt(double(1)+double(.0902)*sqrt(X))-double(1);
        quantity<si::dimensionless, double> a = double(9.06 * 9.06 / 4)
          * pow_hlpr * pow_hlpr / pow(X,b) ;

        quantity<si::dimensionless, double> Av;
        quantity<si::dimensionless, double> Bv;

        if(spherical)
        {
            /// eq 3.1 in @copydetails Khvorostyanov_and_Curry_2002 J. Atmos. Sci 
          Av = a
            * pow(eta_dbl / rhoa_dbl * double(1e4) * si::seconds / si::square_metres, double(1)-double(2)*b)
            * pow(double(4./3) * rho_w<double>() / rhoa_dbl * g<double>() *double(1e2)* si::seconds * si::seconds / si::metres, b) ;
        }
        else
        // nonspherical
        {
          // aspect ratio eq. 3.4
          quantity<si::length, double> lambda_half = double(2.35e-3) * si::metres;
          quantity<si::dimensionless, double> ksi = exp(-r_dbl / lambda_half) 
            + (double(1) - exp(-r_dbl / lambda_half)) / (double(1) + r_dbl / lambda_half);

          // parameters from table 1
          quantity<si::mass_density, double> alfa = 
#if !defined(__NVCC__)
            pi<double>()
#else
            CUDART_PI
#endif
            / double(6) * rho_w<double>() * ksi;    

          // eqs. 2.24, 2.25
          Av = a 
            * pow(eta_dbl / rhoa_dbl * double(1e4) * si::seconds / si::square_metres, double(1)-double(2)*b)
            * pow(double(2.546479) * alfa / rhoa_dbl * g<double>()  * double(1e2) * si::seconds * si::seconds / si::metres  , b);
        }
        Bv = double(3)*b - double(1) ;

        return quantity<si::velocity, real_t>((Av * double(pow(double(2*1e2) * r_dbl/si::metres, Bv)))/double(1e2) * si::metres_per_second);
      }

      // terminal fall velocity at sea level according to Beard (1977)
      // has to be calculated using double prec, on single prec with -use_fast_math it fails on CUDA
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::velocity, real_t> vt_beard77_v0( 
        quantity<si::length, real_t> r //radius
      ) 
      {
#if !defined(__NVCC__)
        using std::pow;
#endif
        // use 3rd degree polynominal for r<20um
        double m_s[4] = {0.105035e2, 0.108750e1, -0.133245, -0.659969e-2};
        // use 7th degree polynominal for r>20um
        double m_l[8] = { 0.65639e1,    -0.10391e1,    -0.14001e1,    -0.82736e0,    -0.34277e0,    -0.83072e-1,    -0.10583e-1,    -0.54208e-3};

        double x = log(2*100*(r / si::metres));
        double y = 0;
        // calc V0 (sea-level velocity)
        if(r <= quantity<si::length, double>(double(20e-6) * si::meters))
          for(int i=0; i<4; ++i)
            y += m_s[i] * pow(x, double(i));
        else
          for(int i=0; i<8; ++i)
            y += m_l[i] * pow(x, double(i));

        return quantity<si::velocity, real_t>((exp(y) / 100.) * si::metres_per_second);
      }

      // terminal fall velocity correction for given altitude according to Beard (1977)
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::dimensionless, real_t> vt_beard77_fact( 
        quantity<si::length, real_t> r, //radius
        quantity<si::pressure, real_t> p, //pressure
        quantity<si::mass_density, real_t> rhoa, //density of air
        quantity<si::dynamic_viscosity, real_t> eta
      ) 
      {
#if !defined(__NVCC__)
        using std::sqrt;
#endif
        using earth::rho_stp; //note: our rho_stp is ca. 1.225 kg/m^3, while Beard uses 1.204
        quantity<si::dynamic_viscosity, real_t> eta_0(1.818e-5 * si::pascals * si::seconds);

        if(r <= quantity<si::length, real_t>(real_t(20e-6) * si::meters))
        {
          using earth::p_stp;
          quantity<si::length, real_t> l_0(6.62e-8 * si::metres);
          quantity<si::length, real_t> l(l_0 * (eta / eta_0) * sqrt(p_stp<real_t>() / p * rho_stp<real_t>() / rhoa));
          return (eta_0 / eta) * (1 + real_t(1.255) * (l / r)) / (1 + real_t(1.255) * (l_0 / r));
        }
        else
        {
          real_t eps_s = (eta_0 / eta) - 1;
          real_t eps_c = sqrt(rho_stp<real_t>() / rhoa) - 1;
          return real_t(1.104) * eps_s + ( (real_t(1.058)*eps_c - real_t(1.104)*eps_s) * (real_t(5.52) + log(2*100 * (r / si::metres)))/real_t(5.01)) +1;
        }
      }
 
      // the exact formula from Beard 1976
      // has to be calculated using double prec, on single prec with -use_fast_math it fails on CUDA
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::velocity, real_t> vt_beard76( 
        quantity<si::length, real_t> r, //radius
        quantity<si::temperature, real_t> T, //temperature
        quantity<si::pressure, real_t> p, //pressure
        quantity<si::mass_density, real_t> rhoa, //density of air
        quantity<si::dynamic_viscosity, real_t> eta
      ) 
      {
#if !defined(__NVCC__)
        using std::pow;
        using std::sqrt;
#endif
        using earth::p_stp;
        using earth::g;
        using moist_air::rho_w;

        if(r <= quantity<si::length, real_t>(real_t(9.5e-6) * si::meters)) //TODO: < 0.5um
        {
          quantity<si::dimensionless, real_t> l = ( real_t(6.62e-8)  * (eta / si::pascals / si::seconds/ real_t(1.818e-5) )  * (p_stp<real_t>() / p)  *  sqrt(real_t(T / si::kelvins) / real_t(293.15)) );
          quantity<si::dimensionless, real_t> C_ac = real_t(1.) + real_t(1.255) * l * si::meters / r;
          return ( (rho_w<real_t>()-rhoa) * g<real_t>() / ( real_t(4.5) * eta) * C_ac * r *r);
        } 

        else if(r <= quantity<si::length, real_t>(real_t(5.035e-4) * si::meters))
        {
          const double b[7] = { -0.318657e1, 0.992696, -0.153193e-2, -0.987059e-3, -0.578878e-3, 0.855176e-4,-0.327815e-5};
          quantity<si::dimensionless, real_t> l = ( real_t(6.62e-8)  * (eta / si::pascals / si::seconds/ real_t(1.818e-5) )  * (p_stp<real_t>() / p)  *  sqrt(real_t(T / si::kelvins) / real_t(293.15)) );
          quantity<si::dimensionless, real_t> C_ac = real_t(1.) + real_t(1.255) * l * si::meters / r;
          quantity<si::dimensionless, real_t> log_N_Da = log( real_t(32./3.) * r * r * r * rhoa * (rho_w<real_t>() - rhoa) * g<real_t>() / eta / eta );
          quantity<si::dimensionless, real_t> Y = 0.;
          for(int i=0; i<7; ++i)
            Y = double(Y) + b[i] * pow(double(log_N_Da), double(i));
          quantity<si::dimensionless, real_t> N_Re = C_ac * exp(double(Y));
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
