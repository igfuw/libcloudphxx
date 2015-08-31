#pragma once

#include <iostream>

#include <blitz/array.h> 

#include <libcloudph++/common/hydrostatic.hpp>
#include <libcloudph++/common/theta_std.hpp>
#include <libcloudph++/common/theta_dry.hpp>
#include <libcloudph++/common/lognormal.hpp>
#include <libcloudph++/common/unary_function.hpp>

#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>

// TODO: relaxation terms still missing

// 8th ICMW case 1 by Wojciech Grabowski)
namespace icmw8_case1
{
  using real_t = float;

  namespace hydrostatic = libcloudphxx::common::hydrostatic;
  namespace theta_std = libcloudphxx::common::theta_std;
  namespace theta_dry = libcloudphxx::common::theta_dry;
  namespace lognormal = libcloudphxx::common::lognormal;

  enum {x, z}; // dimensions

  const quantity<si::temperature, real_t> 
    th_0 = 289 * si::kelvins;
  const quantity<si::dimensionless, real_t> 
    rv_0 = 7.5e-3;
  const quantity<si::pressure, real_t> 
    p_0 = 101500 * si::pascals;
  const quantity<si::velocity, real_t> 
    w_max = real_t(.6) * si::metres_per_second;
  const quantity<si::length, real_t> 
    z_0  = 0    * si::metres,
    Z    = 1500 * si::metres, 
    X    = 1500 * si::metres;
  const quantity<si::time, real_t>
    dt = real_t(1) * si::seconds;

  //aerosol bimodal lognormal dist. 
  const quantity<si::length, real_t>
    mean_rd1 = real_t(.04e-6 / 2) * si::metres,
    mean_rd2 = real_t(.15e-6 / 2) * si::metres;
  const quantity<si::dimensionless, real_t>
    sdev_rd1 = real_t(1.4),
    sdev_rd2 = real_t(1.6);
  const quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t>
    n1_stp = real_t(60e6) / si::cubic_metres,
    n2_stp = real_t(40e6) / si::cubic_metres;

  //aerosol chemical composition parameters (needed for activation)
  // for lgrngn:
  const quantity<si::dimensionless, real_t> kappa = .61; // CCN-derived value from Table 1 in Petters and Kreidenweis 2007
  // for blk_2m:
  const quantity<si::dimensionless, real_t> chem_b = .55; //ammonium sulphate //chem_b = 1.33; // sodium chloride

  //th and rv relaxation time and height
  const quantity<si::time, real_t> tau_rlx = 300 * si::seconds;
  const quantity<si::length, real_t> z_rlx = 200 * si::metres;

  // density profile as a function of altitude
  struct rhod
  {
    real_t operator()(real_t z) const
    {
      quantity<si::pressure, real_t> p = hydrostatic::p(
	z * si::metres, th_0, rv_0, z_0, p_0
      );
      
      quantity<si::mass_density, real_t> rhod = theta_std::rhod(
	p, th_0, rv_0
      );

      return rhod / si::kilograms * si::cubic_metres;
    }

    // to make the rhod() functor accept Blitz arrays as arguments
    BZ_DECLARE_FUNCTOR(rhod);
  };



  /// (similar to eq. 2 in @copydetails Rasinski_et_al_2011, Atmos. Res. 102)
  /// @arg xX = x / X
  /// @arg zZ = z / Z
  real_t psi(real_t xX, real_t zZ) // for computing a numerical derivative
  {
    using namespace boost::math;
    return - sin_pi(zZ) * cos_pi(2 * xX);
  }
  BZ_DECLARE_FUNCTION2_RET(psi, real_t)

  real_t dpsi_dz(real_t xX, real_t zZ)
  {
    using namespace boost::math;
    return - pi<real_t>() / (Z / si::metres) * cos_pi(2 * xX) * cos_pi(zZ);
  }
  BZ_DECLARE_FUNCTION2_RET(dpsi_dz, real_t)

  real_t dpsi_dx(real_t xX, real_t zZ)
  {
    using namespace boost::math;
    return 2 * pi<real_t>() / (X / si::metres) * sin_pi(2 * xX) * sin_pi(zZ);
  }
  BZ_DECLARE_FUNCTION2_RET(dpsi_dx, real_t)

  // function expecting a libmpdata solver parameters struct as argument
  template <class T>
  void setopts(T &params, int nx, int nz)
  {
    params.dt = dt / si::seconds;
    params.dx = (X / si::metres) / (nx-1); 
    params.dz = (Z / si::metres) / (nz-1);
  }

  // function expecting a libmpdata++ solver as argument
  template <class concurr_t>
  void intcond(concurr_t &solver)
  {
    using ix = typename concurr_t::solver_t::ix;

    // helper ondex placeholders
    blitz::firstIndex i;
    blitz::secondIndex j;

    // dx, dy ensuring 1500x1500 domain
    int 
      nx = solver.advectee().extent(x), 
      nz = solver.advectee().extent(z); 
    real_t 
      dx = (X / si::metres) / (nx-1), 
      dz = (Z / si::metres) / (nz-1); 
    real_t A = (w_max / si::metres_per_second) * (nx-1) * dx / pi<real_t>() / real_t(2);

    // constant potential temperature & water vapour mixing ratio profiles
    solver.advectee(ix::th) = (theta_dry::std2dry(th_0, rv_0) / si::kelvins); 
    solver.advectee(ix::rv) = real_t(rv_0);

    // density profile
    solver.g_factor() = rhod()(j * dz);

    // momentum field obtained by numerically differentiating a stream function
    solver.advector(x) = - A * 
    // numerical derivative (see note on div values below)
    (
      psi((i+.5)/(nx-1), (j+.5)/(nz-1))- 
      psi((i+.5)/(nx-1), (j-.5)/(nz-1))  
    ) / dz                       
    // analytical derivative (ditto)
    //dpsi_dz((i+.5)/real_t(nx-1), j/real_t(nz-1))
    * (dt / si::seconds) / dx;  // converting to Courant number

    solver.advector(z) = A * 
    // numerical derivative (max(abs(div)) ~ 5e-10)
    (
      psi((i+.5)/(nx-1), (j+.5)/(nz-1)) - 
      psi((i-.5)/(nx-1), (j+.5)/(nz-1))   
    ) / dx 
    // analytical derivative (max(abs(div)) ~ 3e-5)
    //dpsi_dx(i/real_t(nx-1), (j+.5)/real_t(nz-1))
    * (dt / si::seconds) / dz; // converting to Courant number
  }

  // lognormal aerosol distribution
  template <typename T>
  struct log_dry_radii : public libcloudphxx::common::unary_function<T>
  {
    T funval(const T lnrd) const
    {
      return T((
          lognormal::n_e(mean_rd1, sdev_rd1, n1_stp, quantity<si::dimensionless, real_t>(lnrd)) +
          lognormal::n_e(mean_rd2, sdev_rd2, n2_stp, quantity<si::dimensionless, real_t>(lnrd)) 
        ) * si::cubic_metres
      );
    }

    log_dry_radii *do_clone() const 
    { return new log_dry_radii( *this ); }
  };
};
