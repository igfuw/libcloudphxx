#pragma once

#include "units.hpp"
#include "macros.hpp" 
#include "molar_mass.hpp" 

namespace libcloudphxx
{
  namespace common
  {
    typedef divide_typeof_helper<si::energy, si::temperature>::type energy_over_temperature;
    typedef divide_typeof_helper<si::energy, si::mass>::type energy_over_mass;
    typedef divide_typeof_helper<energy_over_temperature, si::amount>::type energy_over_temperature_amount;
    typedef divide_typeof_helper<energy_over_temperature, si::mass>::type energy_over_temperature_mass;

    typedef multiply_typeof_helper<si::velocity, si::length>::type diffusivity;
    typedef multiply_typeof_helper<si::time, si::area>::type time_area;
    typedef divide_typeof_helper<si::mass, time_area>::type mass_flux;
    typedef multiply_typeof_helper<energy_over_mass, mass_flux>::type energy_flux;
    typedef divide_typeof_helper<si::temperature, si::length>::type temperature_gradient;
    typedef divide_typeof_helper<energy_flux, temperature_gradient>::type thermal_conductivity;

    namespace moist_air
    {
      // specific heat capacities
      libcloudphxx_const(energy_over_temperature_mass, c_pd, 1005, si::joules / si::kilograms / si::kelvins) // dry air
      libcloudphxx_const(energy_over_temperature_mass, c_pv, 1850, si::joules / si::kilograms / si::kelvins) // water vapour
      libcloudphxx_const(energy_over_temperature_mass, c_pw, 4218, si::joules / si::kilograms / si::kelvins) // liquid water
      libcloudphxx_const(energy_over_temperature_mass, c_pi, 2114, si::joules / si::kilograms / si::kelvins) // ice // ice

      // molar masses
      libcloudphxx_const(mass_over_amount, M_d, 0.02897, si::kilograms / si::moles) // dry air (Curry & Webster / Seinfeld & Pandis)
      libcloudphxx_const_derived(mass_over_amount, M_v, molar_mass::M_H<real_t>() + molar_mass::M_OH<real_t>()) // water vapour
      libcloudphxx_const_derived(si::dimensionless, eps, M_v<real_t>() / M_d<real_t>()) // aka epsilon

      // universal gas constant (i.e. the Boltzmann times the Avogadro constants)
      // source: http://physics.nist.gov/cgi-bin/cuu/Value?r
      //         http://glossary.ametsoc.org/wiki/Gas_constant
      //         Mohr, P. J., B. N. Taylor, and D. B. Newell, 2012: CODATA recommended values of the fundamental physical constants: 2010. J. Phys. Chem. Ref. Data, 41, 043109, doi:10.1063/1.4724320. 
      libcloudphxx_const(energy_over_temperature_amount, kaBoNA, 8.3144621, si::joules / si::kelvin / si::mole)

      // gas constants
      libcloudphxx_const_derived(energy_over_temperature_mass, R_d, kaBoNA<real_t>() / M_d<real_t>()) // dry air
      libcloudphxx_const_derived(energy_over_temperature_mass, R_v, kaBoNA<real_t>() / M_v<real_t>()) // water vapour

      // Exner function exponent for dry air
      libcloudphxx_const_derived(si::dimensionless, R_d_over_c_pd, R_d<real_t>() / c_pd<real_t>())

      // water density
      libcloudphxx_const(si::mass_density, rho_w, 1e3, si::kilograms / si::cubic_metres)
      // ice density
      libcloudphxx_const(si::mass_density, rho_i, 910, si::kilograms / si::cubic_metres)
      // graupel density used for ice B (from Grabowski 1999)
      libcloudphxx_const(si::mass_density, rho_ib, 400, si::kilograms / si::cubic_metres)

      // mixing rule for extensive quantitites (i.e. using mass mixing ratio)
      template <typename real_t, typename quant>
      quant constexpr mix(
	const quant &dry, 
	const quant &vap, 
	const quantity<si::dimensionless, real_t> &r
      ) {
	return (dry + r * vap) / (1 + r);
      }

      // gas constant for moist air
      template <typename real_t>
      quantity<energy_over_temperature_mass, real_t> R(
	const quantity<si::dimensionless, real_t> &r
      ) {
	return mix(R_d<real_t>(), R_v<real_t>(), r);
      }
     
      // specific heat capacity for moist air
      template <typename real_t>
      quantity<energy_over_temperature_mass, real_t> c_p(
	const quantity<si::dimensionless, real_t> &r
      ) {
	return mix(c_pd<real_t>(), c_pv<real_t>(), r);
      }
      
      // water vapour partial pressure as a function of mixing ratio
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::pressure, real_t> p_v(
	const quantity<si::pressure, real_t> &p,
	const quantity<si::dimensionless, real_t> &r
      ) {
	return p * r / (r + eps<real_t>());
      }

      //vapour diffusivity in air (see Properties of air, Tracy, Welch & Porter 1980)
      libcloudphxx_const(diffusivity, D_0, 2.26e-5, si::square_metres / si::seconds) 

      template<typename real_t>
      BOOST_GPU_ENABLED
      quantity<diffusivity, real_t> D(
        const quantity<si::temperature, real_t> &T, 
        const quantity<si::pressure, real_t> &p
      ) { 
        const quantity<si::pressure, real_t> 
          p_0 = real_t(100000) * si::pascal;
        const quantity<si::temperature, real_t> 
          T_0 = real_t(273.15) * si::kelvin;

#if !defined(__NVCC__)
        using std::pow;
#endif
 
        return D_0<real_t>() * pow(T / T_0, real_t(1.81)) * (p_0 / p); 
      }   

      // thermal conductivity of air
      libcloudphxx_const(thermal_conductivity, K_0, 2.4e-2, si::joules / si::metres / si::seconds / si::kelvins) 
    };
  };
};
