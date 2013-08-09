#pragma once

#include <libcloudph++/common/units.hpp>
#include <libcloudph++/common/macros.hpp>
#include <libcloudph++/common/const_cp.hpp>

namespace libcloudphxx
{
  namespace common
  {
    namespace maxwell_mason
    {
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<divide_typeof_helper<si::area, si::time>::type, real_t> rdrdt(
	const quantity<si::mass_density, real_t> rho_v, // ambient water vapour density
	const quantity<si::temperature, real_t> T, // ambient temperature
	const quantity<si::pressure, real_t> p, // ambient pressure
	const quantity<si::dimensionless, real_t> RH, // p_v/p_vs = relative humidity
	const quantity<si::dimensionless, real_t> a_w, // water activity
	const quantity<si::dimensionless, real_t> klvntrm // the Kelvin term
      )
      {
        using namespace moist_air;

	quantity<divide_typeof_helper<si::energy, si::mass>::type, real_t> l_v = const_cp::l_v<real_t>(T);
	return 
	  (real_t(1) - a_w * klvntrm / RH)
	  / rho_w<real_t>() 
	  / ( 
	    real_t(1) / D<real_t>(T, p) / rho_v
	    +
	    l_v / K_0<real_t>() / RH / T * (l_v / R_v<real_t>() / T - real_t(1))
	  )   
	;   
      }
    };
  };
};
