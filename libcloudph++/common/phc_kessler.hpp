#pragma once
#include <libcloudph++/common/phc.hpp>

namespace phc
{
// TODO: namespace kessler?

  //Kessler autoconversion
  //eq. 5a in Grabowski & Smolarkiewicz 1996
  phc_declare_const_macro(k_autoconv_default, .001, si::hertz)
  phc_declare_const_macro(rc_thresh_default, .0005, si::dimensionless()) 

  phc_declare_funct_macro quantity<divide_typeof_helper<mixing_ratio, si::time>::type, real_t> autoconversion_rate(
    quantity<mixing_ratio,  real_t> rc,
    quantity<mixing_ratio,  real_t> rc_thresh = rc_thresh_default<real_t>(),
    quantity<si::frequency, real_t> k_autoconv = k_autoconv_default<real_t>()
  )
  {
    return k_autoconv * std::max( real_t(0) * si::dimensionless(), rc - rc_thresh);
  }

  //Kessler collection
  //eq. 5b in Grabowski & Smolarkiewicz 1996
  phc_declare_const_macro(k_2, 2.2, si::hertz)

  phc_declare_funct_macro quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> collection_rate(
    quantity<mixing_ratio,  real_t> rc,
    quantity<mixing_ratio,  real_t> rr
  )
  {
    return k_2<real_t>() * rc * pow(rr, .875);
  }

  // Kessler evaporation rate
  // eq. 5c in Grabowski & Smolarkiewicz 1996 (multiplied by rho!)
  phc_declare_funct_macro quantity<divide_typeof_helper<si::mass_density, si::time>::type, real_t> evaporation_rate( 
    quantity<mixing_ratio, real_t> rv,
    quantity<mixing_ratio, real_t> rvs,
    quantity<si::mass_density, real_t> rhod_rr,
    quantity<si::pressure, real_t> p
  )
  {
    return (1 - rv / rvs) 
      * (1.6 + 124.9 * pow(1e-3 * (rhod_rr * si::cubic_metres / si::kilograms), .2046)) // ventilation factor
      * pow(1e-3 * (rhod_rr * si::cubic_metres / si::kilograms), .525) 
      / (5.4e2 + 2.55e5 * (1. / (p / si::pascals) / rvs)) 
      / si::seconds * si::kilograms / si::cubic_metres;
  }

  // Kessler/Beard terminal velocity
  // eq. 5d in Grabowski & Smolarkiewicz 1996	
  phc_declare_const_macro(vterm_A, 36.34, si::metre_per_second)
  phc_declare_const_macro(vterm_B, 1e-3, si::cubic_metres / si::kilograms)

  phc_declare_funct_macro quantity<si::velocity, real_t> v_term(
    quantity<si::mass_density, real_t> rho_r,
    quantity<si::mass_density, real_t> rho_d,
    quantity<si::mass_density, real_t> rho_d0
  )
  {
    return vterm_A<real_t>() * real_t(pow(rho_r * vterm_B<real_t>(), .1346) * sqrt(rho_d0 / rho_d));
  }

};
