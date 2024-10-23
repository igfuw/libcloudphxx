/** @file
  * @copyright University of Warsaw
  * @brief single-moment bulk parameterisation formulae (Kessler)
  *   from @copybrief bib::Grabowski_and_Smolarkiewicz_1996
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include "../common/moist_air.hpp"

namespace libcloudphxx
{
  namespace blk_1m
  {
    namespace formulae
    {
      //Kessler autoconversion
      //eq. 5a in Grabowski & Smolarkiewicz 1996
      //(the default k_autoconv is 0.001 1/s)
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> autoconversion_rate(
        const quantity<si::dimensionless, real_t> &rc,
        const quantity<si::dimensionless, real_t> rc_thresh,
        const quantity<si::frequency, real_t> k_autoconv
      ) {
        return k_autoconv * std::max( real_t(0) * si::dimensionless(), rc - rc_thresh);
      }

      //Kessler collection
      //eq. 5b in Grabowski & Smolarkiewicz 1996
      libcloudphxx_const(si::frequency, k_2, 2.2, si::hertz)

      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> collection_rate(
        const quantity<si::dimensionless, real_t> &rc,
        const quantity<si::dimensionless, real_t> &rr
      ) {
        return k_2<real_t>() * rc * std::pow(rr, real_t(.875));
      }

      // Kessler evaporation rate
      // eq. 5c in Grabowski & Smolarkiewicz 1996 (multiplied by rho!)
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> evaporation_rate( 
        quantity<si::dimensionless, real_t> rv,
        quantity<si::dimensionless, real_t> rvs,
        quantity<si::dimensionless, real_t> rr,
        quantity<si::mass_density, real_t> rhod,
        quantity<si::pressure, real_t> p
      )
      {
        return 
          (1 - rv / rvs) / rhod
          * (
            real_t(1.6) 
            + real_t(124.9) * std::pow( 
              real_t(1e-3) * rhod * rr * si::cubic_metres / si::kilograms,
              real_t(.2046)
            ) 
          ) // ventilation factor TODO- move to ventil.hpp
          * std::pow(
            real_t(1e-3) * rhod * rr * si::cubic_metres / si::kilograms, 
            real_t(.525)
          ) 
          / (real_t(5.4e2) 
          + real_t(2.55e5) 
          * (real_t(1) / (p / si::pascals) / rvs)) 
          / si::seconds * si::kilograms / si::cubic_metres;
      }

      // Kessler/Beard terminal velocity
      // eq. 5d in Grabowski & Smolarkiewicz 1996        
      libcloudphxx_const(si::velocity, vterm_A, 36.34, si::metre_per_second)

      using inverse_density = divide_typeof_helper<si::dimensionless,si::mass_density>::type;
      libcloudphxx_const(inverse_density, vterm_B, 1e-3, si::cubic_metres / si::kilograms)

      template <typename real_t>
      quantity<si::velocity, real_t> v_term(
        const quantity<si::dimensionless, real_t> &rr,
        const quantity<si::mass_density, real_t> &rhod,
        const quantity<si::mass_density, real_t> &rhod_0
      ) {
        return 
          vterm_A<real_t>() 
          * real_t(std::pow(
            (rhod * rr * vterm_B<real_t>()), 
            real_t(.1346)
          ) 
          * sqrt(rhod_0 / rhod)
        );
      }

      //ice A homogeneous nucleation 1
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> hom_A_nucleation_1(
        const quantity<si::dimensionless, real_t> &rv,
        const quantity<si::dimensionless, real_t> &rvs,
        const quantity<si::dimensionless, real_t> &rvsi,
        const quantity<si::kelvins, real_t> &th,
        const quantity<si::pascals, real_t> &p,
        const quantity<si::time, real_t> &dt
      ) {
        const quantity<si::kelvin, real_t> T = std::pow((p/real_t(100000)/si::pascal), 0.286);
        const quantity<si::time, real_t> taunuc = dt; // time scale for nucleation
        //homogeneous part
        if (T < real_t(233.16) * si::kelvins) { // only if T < -40 C
        real_t beta = (T > real_t(213.16)*si::kelvins) ? real_t(0.1)+real_t(0.9)*(T-real_t(213.16)*si::kelvins)/(real_t(20)*si::kelvins) : real_t(0.1);
        quantity<si::dimensionless, real_t> rv_adj = beta * rvs + (1-beta) * rvsi;
        //nucleation rate:
        quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> hom1 = (real_t(1)-std::exp(-dt/taunuc)) * std::max(real_t(0), rv-rv_adj) / dt;
        }
        else {
          quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> hom1 = real_t(0);
        }
        return hom1;
      }

      //ice A homogeneous nucleation 2
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> hom_A_nucleation_2(
        const quantity<si::dimensionless, real_t> &rc,
        const quantity<si::kelvins, real_t> &th,
        const quantity<si::pascals, real_t> &p,
        const quantity<si::time, real_t> &dt
      ) {
        const quantity<si::kelvin, real_t> T = std::pow((p/real_t(100000)/si::pascal), 0.286);
        const quantity<si::time, real_t> taunuc = dt; // time scale for nucleation
        //homogeneous part
        if (T < real_t(233.16) * si::kelvins) { // only if T < -40 C
        //nucleation rate:
        quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> hom2 = (real_t(1)-std::exp(-dt/taunuc)) * rc / dt;
        }
        else {
          quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> hom2 = real_t(0);
        }
        return hom2;
      }

      //ice A heterogeneous nucleation
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> het_A_nucleation(
        const quantity<si::dimensionless, real_t> &ra,
        const quantity<si::dimensionless, real_t> &rc,
        const quantity<si::kelvins, real_t> &th,
        const quantity<si::pascals, real_t> &p,
        const quantity<si::mass_density, real_t> &rhod_0,
        const quantity<si::time, real_t> &dt
      ) {
        const quantity<si::kelvin, real_t> T = std::pow((p/real_t(100000)/si::pascal), 0.286);
        const quantity<si::time, real_t> taunuc = dt; // time scale for nucleation
        //heterogeneous part
        //calculating the mass of small ice A particle:
        real_t IWCS = std::min(real_t(1e-3), rhod_0*ra/si::mass_density, real_t(2.52e-4)*std::pow(real_t(1e3)*rhod_0*ra/si::mass_density,real_t(0.837)));
        real_t rho_i = real_t(916.8); //ice density
        real_t alpha = std::max(real_t(2.2e4), real_t(-4.99e3)-real_t(4.94e4)*std::log10(real_t(1e3)*IWCS));
        quantity<si::mass, real_t> m_as = real_t(2*std::numbers::pi)*rho_i/std::pow(alpha,3) *si::mass;
        //calculating the mass of large ice A particle:
        real_t IWCL = rhod_0*ra - IWCS;
        real_t a_mu = real_t(5.2) + real_t(1.3e-3)*(T/si::kelvin-273.16);
        real_t b_mu = real_t(0.026) - real_t(1.2e-3)*(T/si::kelvin-273.16);
        real_t a_sigma = real_t(0.47) + real_t(2.1e-3)*(T/si::kelvin-273.16);
        real_t b_sigma = real_t(0.018) - real_t(2.1e-4)*(T/si::kelvin-273.16);
        real_t mu = a_mu + b_mu*std::log10(real_t(1e3)*IWCL);
        real_t sigma = a_sigma + b_sigma*std::log10(real_t(1e3)*IWCL);
        quantity<si::mass, real_t> m_al = real_t(1.67e17*pi<real_t>())*rho_i*std::exp(3*mu +9/2*std::pow(sigma,2))*si::mass;
        //mass of average ice A particle
        real_t delta = IWCS/(ra*rhod_0/si::mass_density);
        quantity<si::mass, real_t> ma = delta*m_as + (real_t(1)-delta)*m_al;
        //concentration of ice nuclei:
        quantity<divide_typeof_helper<si::dimensionless, si::cubic_meters>::type, real_t> N_in = std::min(real_t(1e5), real_t(1e-2)*std::exp(real_t(0.6)*(real_t(273.16)-T/si::kelvins)));
        //nucleation rate:
        quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> het = (real_t(1)-std::exp(-dt/taunuc)) * std::min(rc, std::max(real_t(0), N_in*std::max(real_t(1e-12)*si::mass, ma)/rhod_0 - ra)) / dt;
        return het;
      }

    };
  };
};
