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
};
