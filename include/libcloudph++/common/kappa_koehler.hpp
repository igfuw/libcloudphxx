#pragma once

#include "units.hpp"
#include "macros.hpp"
#include "kelvin_term.hpp"
#include "detail/toms748.hpp"

namespace libcloudphxx
{
  namespace common
  {
    namespace kappa_koehler
    {
      /// @brief equilibrium wet radius to the third power, with the Kelvin term discarded, for a given:
      /// @arg dry radius to the third power
      /// @arg the solubility parameter kappa
      /// @arg ratio of abmient vapour density/pressure to saturation vapour density/pressure for pure water
      /// 
      /// the formula stems from applying the kappa-Koehler relation 
      /// (eq. 6 in @copydetails Petters_and_Kreidenweis_2007) to a stationarity
      /// condition for a vapour diffusion equation, which translates to
      /// zero differece of vapour density: rho_ambient - rho_surface = 0
      /// 
      /// since rho_surface = rho_surface_pure_water * a(r_w, r_d, kappa)
      /// one can derive r_w as a function of r_d, kappa and the ratio
      /// of abmient and surface vapour densities
      /// 
      /// for the kappa-Koehler parameterisation rw3 is linear with rd3
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::volume, real_t> rw3_eq_nokelvin(
	quantity<si::volume, real_t> rd3, 
	quantity<si::dimensionless, real_t> kappa,
	quantity<si::dimensionless, real_t> RH
      )   
      {   
        assert(RH > .05); // kappa-koehler assumes well dissolved matter
        assert(RH < 1.); // no equilibrium over RH=100%
        assert(kappa > 0); // pure-water case left out
	return rd3 * (1 - RH * (1 - kappa)) / (1 - RH);
      }   

      /// @brief activity of water in solution (eqs. 1,6) in @copydetails Petters_and_Kreidenweis_2007
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::dimensionless, real_t> a_w(
	quantity<si::volume, real_t> rw3,
	quantity<si::volume, real_t> rd3,
	quantity<si::dimensionless, real_t> kappa
      )
      {
        assert(kappa > 0); // pure-water case left out
	return (rw3 - rd3) / (rw3 - rd3 * (real_t(1) - kappa));
      }

      // functor to be passed to the minimisation functions used below
      namespace detail
      {
        // equilibrium radius
        template <typename real_t>
	struct rw3_eq_minfun 
	{   
	  const quantity<si::dimensionless, real_t> RH;
	  const quantity<si::volume, real_t> rd3;
	  const quantity<si::dimensionless, real_t> kappa;
	  const quantity<si::temperature, real_t> T;

          // ctor
          BOOST_GPU_ENABLED
          rw3_eq_minfun(
	    const quantity<si::dimensionless, real_t> &RH,
	    const quantity<si::volume, real_t> &rd3,
	    const quantity<si::dimensionless, real_t> &kappa,
	    const quantity<si::temperature, real_t> &T
          ) : RH(RH), rd3(rd3), kappa(kappa), T(T) {}
      
          BOOST_GPU_ENABLED
	  real_t operator()(real_t rw3) const
	  {
#if !defined(__NVCC__)
	    //using std::pow;
#endif
	    return this->RH
	      - a_w(rw3 * si::cubic_metres, this->rd3, this->kappa) 
	      * kelvin::klvntrm(cbrt(rw3) * si::metres, this->T); 
	  }
	};  

        // critical radius
        template <typename real_t>
        struct rw3_cr_minfun
        {
	  const quantity<si::volume, real_t> rd3;
	  const quantity<si::dimensionless, real_t> kappa;
	  const quantity<si::temperature, real_t> T;

          // ctor
          BOOST_GPU_ENABLED
          rw3_cr_minfun(
	    const quantity<si::volume, real_t> &rd3,
	    const quantity<si::dimensionless, real_t> &kappa,
	    const quantity<si::temperature, real_t> &T
          ) : rd3(rd3), kappa(kappa), T(T) {}

          BOOST_GPU_ENABLED
          real_t operator()(real_t rw3) const
          {
#if !defined(__NVCC__)
	    //using std::pow;
#endif
            return (kelvin::A(T) 
              * (rd3 - rw3 * si::cubic_metres) 
              * ((kappa - 1) * rd3 + rw3 * si::cubic_metres) 
              + 3 * kappa * rd3 * rw3 * cbrt(rw3) * si::cubic_metres * si::metres
            ) / si::cubic_metres / si::cubic_metres / si::metres;
          }
        };
      };


      // @brief equilibrium wet radius to the third power for a given:
      /// @arg dry radius to the third power
      /// @arg the solubility parameter kappa
      /// @arg ratio of abmient vapour density/pressure to saturation vapour density/pressure for pure water
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::volume, real_t> rw3_eq(
	quantity<si::volume, real_t> rd3, 
	quantity<si::dimensionless, real_t> kappa,
	quantity<si::dimensionless, real_t> RH,
	quantity<si::temperature, real_t> T
      )   
      {   
        assert(RH < 1); // no equilibrium over RH=100%
        assert(kappa > 0); // pure-water case left out

        return common::detail::toms748_solve(
	  detail::rw3_eq_minfun<real_t>(RH, rd3, kappa, T), // the above-defined functor
	  real_t(rd3 / si::cubic_metres), // min
	  real_t(rw3_eq_nokelvin(rd3, kappa, RH) / si::cubic_metres) // max
	) * si::cubic_metres;
      }

      // critical radius, i.e. radius for which the Koehler curve reaches maximum
      // formula obtained by analytically derivating the kappa-Koehler curve:
      // enter "derivative of (x^3-r^3)/(x^3-r^3*(1-k))*exp(A/x)" at wolframalpha.com 
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::volume, real_t> rw3_cr(
	quantity<si::volume, real_t> rd3, 
	quantity<si::dimensionless, real_t> kappa,
	quantity<si::temperature, real_t> T
      )   
      {   
        assert(kappa > 0); // pure-water case left out

        quantity<si::volume, double> rd3_dbl = static_cast<quantity<si::volume, double> >(rd3);
        quantity<si::temperature, double> T_dbl = static_cast<quantity<si::temperature, double> >(T);

        return real_t(common::detail::toms748_solve(
	  detail::rw3_cr_minfun<double>(rd3_dbl, double(kappa), T_dbl), // the above-defined functor
	  double(1e0 * (rd3 / si::cubic_metres)), // min
	  double(1e8 * (rd3 / si::cubic_metres))  // max
	)) * si::cubic_metres;
      }

      // critical supersaturation, i.e. S(r_cr)
      template <typename real_t>
      BOOST_GPU_ENABLED
      quantity<si::dimensionless, real_t> S_cr(
        quantity<si::volume, real_t> rd3,
        quantity<si::dimensionless, real_t> kappa,
        quantity<si::temperature, real_t> T
      )
      {
        assert(kappa > 0); // pure-water case left out
#if !defined(__NVCC__)
        //using std::pow;
#endif
        quantity<si::volume, real_t> rw3 = rw3_cr(rd3, kappa, T);

        return a_w(rw3, rd3, kappa) * kelvin::klvntrm(
          cbrt(rw3 / si::cubic_metres) * si::metres, 
          T
        );
      }
    };
  };
};
