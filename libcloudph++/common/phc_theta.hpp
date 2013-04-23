#pragma once

#include <libcloudph++/common/phc.hpp>
#include <libcloudph++/common/phc_const_cp.hpp>

namespace phc
{
  // Exner function exponent for moist air
  phc_declare_funct_macro quantity<si::dimensionless, real_t> R_over_c_p(
    quantity<mixing_ratio, real_t> r
  )
  {
    return R<real_t>(r) / c_p<real_t>(r);
  }


  // Exner function for dry air
  phc_declare_funct_macro quantity<si::dimensionless, real_t> exner(
    quantity<si::pressure, real_t> p
  )
  {
    return pow(p / p_1000<real_t>(), R_d_over_c_pd<real_t>());
  }


  // Exner function for moist air
  phc_declare_funct_macro quantity<si::dimensionless, real_t> exner(
    quantity<si::pressure, real_t> p,
    quantity<mixing_ratio, real_t> r
  )
  {
    return pow(p / p_1000<real_t>(), R_over_c_p<real_t>(r));
  }

  // dry air density as a function of p, theta and rv
  phc_declare_funct_macro quantity<si::mass_density, real_t> rhod(
    quantity<si::pressure, real_t> p,
    quantity<si::temperature, real_t> th,
    quantity<phc::mixing_ratio, real_t> rv
  )
  {
    return (p - phc::p_v<real_t>(p, rv)) /
      (phc::exner<real_t>(p, rv) * phc::R_d<real_t>() * th);
  }

  // temperature as a function theta, pressure and water vapour mixing ratio for moist air
  phc_declare_funct_macro quantity<si::temperature, real_t> T(
    quantity<si::temperature, real_t> th, // TODO: const
    quantity<si::pressure, real_t> p, // TODO: const
    quantity<mixing_ratio, real_t> r // TODO: const
  )
  {
    return th * exner<real_t>(p, r);
  }


  // pressure as a function of "theta times dry air density" and water vapour mixing ratio
  phc_declare_funct_macro quantity<si::pressure, real_t> p(
    const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> rhod_th,
    quantity<mixing_ratio, real_t> r // TODO: const
  )
  {
    return p_1000<real_t>() * real_t(pow(
      (rhod_th * R_d<real_t>())
        / p_1000<real_t>() * (real_t(1) + r / eps<real_t>()),
      1 / (1 - R_over_c_p(r))
    ));
  }

  // dtheta^star_drv from First Law for theta^star
  // TODO: document the derivation!
  phc_declare_funct_macro quantity<si::temperature, real_t> dtheta_drv(
    const quantity<si::temperature, real_t> T,
    const quantity<si::pressure, real_t> p,
    const quantity<mixing_ratio, real_t> r,
    const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t> rhod_th,
    const quantity<si::mass_density, real_t> rhod
  )
  {
    return - rhod_th / rhod * (
      // the 'liquid water' term
      l_v<real_t>(T)
        / real_t(pow(1 + r, 2))
        / c_p(r)
        / T
/*
+
// the 'virtual' term (TODO: as an option!)
log(p / p_1000<real_t>())
* R_d_over_c_pd<real_t>()
* (real_t(1) / eps<real_t>() - real_t(1) / ups<real_t>())
* real_t(pow(real_t(1) + r / ups<real_t>(),-2))
*/
    );
  }
};
