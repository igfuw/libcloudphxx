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

// setup taken from 8th ICMW case 1 by Wojciech Grabowski
namespace config
{
  using real_t = float;

  namespace hydrostatic = libcloudphxx::common::hydrostatic;
  namespace theta_std = libcloudphxx::common::theta_std;
  namespace theta_dry = libcloudphxx::common::theta_dry;
  namespace lognormal = libcloudphxx::common::lognormal;

  enum {x, z}; // dimensions

  class setup_t
  {
    public:
    quantity<si::temperature, real_t> th_0;
    quantity<si::dimensionless, real_t> rv_0;
    quantity<si::pressure, real_t> p_0;
    quantity<si::velocity, real_t> w_max;
    quantity<si::length, real_t>  z_0, Z, X;
    quantity<si::time, real_t> dt;

    //aerosol bimodal lognormal dist. 
    quantity<si::length, real_t> mean_rd1, mean_rd2;
    quantity<si::dimensionless, real_t> sdev_rd1, sdev_rd2;
    quantity<power_typeof_helper<si::length, static_rational<-3>>::type, real_t> n1_stp, n2_stp;

    //aerosol chemical composition parameters (needed for activation)
    // for lgrngn:
    quantity<si::dimensionless, real_t> kappa; // CCN-derived value from Table 1 in Petters and Kreidenweis 2007
    // for blk_2m:
    quantity<si::dimensionless, real_t> chem_b; //ammonium sulphate //chem_b = 1.33; // sodium chloride

    //th and rv relaxation time and height
    quantity<si::time, real_t> tau_rlx;
    quantity<si::length, real_t> z_rlx;
  };

  // lognormal aerosol distribution
  template <typename T>
  struct log_dry_radii : public libcloudphxx::common::unary_function<T>
  {
    setup_t setup;
    log_dry_radii(const setup_t &setup):
      setup(setup) {}

    T funval(const T lnrd) const
    {
      return T((
          lognormal::n_e(setup.mean_rd1, setup.sdev_rd1, setup.n1_stp, quantity<si::dimensionless, real_t>(lnrd)) +
          lognormal::n_e(setup.mean_rd2, setup.sdev_rd2, setup.n2_stp, quantity<si::dimensionless, real_t>(lnrd)) 
        ) * si::cubic_metres
      );
    }

    log_dry_radii *do_clone() const 
    { return new log_dry_radii( *this ); }
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
/*
  struct dpsi_dz
  {
    setup_t setup;
    dpsi_dz(const setup_t &setup):
      setup(setup) {}

    real_t operator()(real_t xX, real_t zZ) // for computing a numerical derivative
    {
      using namespace boost::math;
      return - pi<real_t>() / (setup.Z / si::metres) * cos_pi(2 * xX) * cos_pi(zZ);
    }
    BZ_DECLARE_FUNCTOR(dpsi_dz)
  }
  struct dpsi_dx
  {
    setup_t setup;
    dpsi_dx(const setup_t &setup):
      setup(setup) {}

    real_t operator()(real_t xX, real_t zZ) // for computing a numerical derivative
    {
      using namespace boost::math;
      return 2 * pi<real_t>() / (setup.X / si::metres) * sin_pi(2 * xX) * sin_pi(zZ);
    }
    BZ_DECLARE_FUNCTOR(dpsi_dx)
  }
*/
  // density profile as a function of altitude
  struct rhod
  {
    setup_t setup;
    rhod(const setup_t &setup):
      setup(setup)
      {}

    real_t operator()(real_t z) const
    {
      quantity<si::pressure, real_t> p = hydrostatic::p(
	z * si::metres, setup.th_0, setup.rv_0, setup.z_0, setup.p_0
      );
      
      quantity<si::mass_density, real_t> rhod = theta_std::rhod(
	p, setup.th_0, setup.rv_0
      );

      return rhod / si::kilograms * si::cubic_metres;
    }

    // to make the rhod() functor accept Blitz arrays as arguments
    BZ_DECLARE_FUNCTOR(rhod);
  };

  // function expecting a libmpdata solver parameters struct as argument
  template <class T>
  void setopts(T &params, int nx, int nz, setup_t &setup)
  {
    params.dt = setup.dt / si::seconds;
    params.dx = (setup.X / si::metres) / (nx-1); 
    params.dz = (setup.Z / si::metres) / (nz-1);
    params.setup = setup;
  }

  // function expecting a libmpdata++ solver as argument
  template <class concurr_t>
  void intcond(concurr_t &solver, const setup_t &setup)
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
      dx = (setup.X / si::metres) / (nx-1), 
      dz = (setup.Z / si::metres) / (nz-1); 
    real_t A = (setup.w_max / si::metres_per_second) * (nx-1) * dx / pi<real_t>() / real_t(2);

    // constant potential temperature & water vapour mixing ratio profiles
    solver.advectee(ix::th) = (theta_dry::std2dry(setup.th_0, setup.rv_0) / si::kelvins); 
    solver.advectee(ix::rv) = real_t(setup.rv_0);

    // density profile
    solver.g_factor() = rhod(setup)(j * dz);

    // momentum field obtained by numerically differentiating a stream function
    solver.advector(x) = - A * 
    // numerical derivative (see note on div values below)
    (
      psi((i+.5)/(nx-1), (j+.5)/(nz-1))- 
      psi((i+.5)/(nx-1), (j-.5)/(nz-1))  
    ) / dz                       
    // analytical derivative (ditto)
    //dpsi_dz((i+.5)/real_t(nx-1), j/real_t(nz-1))
    * (setup.dt / si::seconds) / dx;  // converting to Courant number

    solver.advector(z) = A * 
    // numerical derivative (max(abs(div)) ~ 5e-10)
    (
      psi((i+.5)/(nx-1), (j+.5)/(nz-1)) - 
      psi((i-.5)/(nx-1), (j+.5)/(nz-1))   
    ) / dx 
    // analytical derivative (max(abs(div)) ~ 3e-5)
    //dpsi_dx(i/real_t(nx-1), (j+.5)/real_t(nz-1))
    * (setup.dt / si::seconds) / dz; // converting to Courant number
  }
};
