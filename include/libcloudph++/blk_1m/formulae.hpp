/** @file
  * @copyright University of Warsaw
  * @brief single-moment bulk parameterisation formulae (Kessler)
  *   from @copybrief bib::Grabowski_and_Smolarkiewicz_1996, Grabowski 1999
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include "../common/moist_air.hpp"
#include "../common/vterm.hpp"

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
        return k_autoconv * std::max(real_t(0) * si::dimensionless(), rc - rc_thresh);
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

      // Kessler rain evaporation rate
      // eq. 5c in Grabowski & Smolarkiewicz 1996 (multiplied by rho!)
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> evaporation_rate(
        quantity<si::dimensionless, real_t> rv,
        quantity<si::dimensionless, real_t> rvs,
        quantity<si::dimensionless, real_t> rr,
        quantity<si::mass_density, real_t> rhod,
        quantity<si::pressure, real_t> p
      ) {
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

      // Kessler/Beard rain terminal velocity
      // eq. 5d in Grabowski & Smolarkiewicz 1996
      libcloudphxx_const(si::velocity, vterm_A, 36.34, si::metre_per_second)

      using inverse_density = divide_typeof_helper<si::dimensionless, si::mass_density>::type;
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

      // slope of Marshall Palmer distribution for rain
      // eq. A.1 in Grabowski 1999
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::length>::type, real_t> lambda_rain(
        const quantity<si::dimensionless, real_t> &rr,
        const quantity<si::mass_density, real_t> &rhod_0
      ) {
        using namespace common;
        real_t N_0r = real_t(1e7);
        return std::pow(pi<real_t>() * moist_air::rho_w<real_t>() * N_0r / (rhod_0 * rr), real_t(0.25)) / si::meters;
      }

      // mean mass of ice A particle
      // eq. A.7 - A.15a in Grabowski 1999
      template <typename real_t>
      quantity<si::mass, real_t> mass_a(
        const quantity<si::dimensionless, real_t> &ria,
        const quantity<si::temperature, real_t> &T,
        const quantity<si::mass_density, real_t> &rhod_0
      ) {
        using namespace common;
        real_t tempc = T / si::kelvin - real_t(273.16);
        real_t IWC = std::max(real_t(1.e-9) * si::dimensionless(), rhod_0 / si::kilograms * si::cubic_meters * ria);
        real_t IWCS = std::max(real_t(1.e-9),
                               std::min(std::min(real_t(1.e-3), IWC),
                                        real_t(2.52e-4) * std::pow((IWC / real_t(1.e-3)), real_t(0.837))));
        real_t IWCL = std::max(real_t(1.e-9), IWC - IWCS);
        //mass of small ice A:
        real_t alpha = std::max(real_t(1.e3), real_t(4.99e3) - real_t(4.94e4) * std::log10(IWCS / real_t(1.e-3)));
        real_t m_as = real_t(6.28) * moist_air::rho_i<real_t>() / si::kilograms * si::cubic_meters / std::pow(
          alpha, real_t(3));
        //mass of large ice A:
        real_t ami = real_t(5.20) + real_t(1.3e-3) * tempc;
        real_t bmi = real_t(0.026) - real_t(1.2e-3) * tempc;
        real_t asi = real_t(0.47) + real_t(2.1e-3) * tempc;
        real_t bsi = real_t(0.018) - real_t(2.1e-4) * tempc;
        real_t alorat = std::log10(IWCL / real_t(1.e-3));
        real_t miu = std::max(real_t(4.6), std::min(real_t(5.4), ami + bmi * alorat));
        real_t sig = std::max(real_t(0.), std::min(real_t(0.5), asi + bsi * alorat));
        real_t exp_mas = real_t(3.) * miu + real_t(4.5) * std::pow(sig, 2);
        real_t m_al = real_t(5.24e-19) * moist_air::rho_i<real_t>() / si::kilograms * si::cubic_meters *
          std::exp(exp_mas);
        //mass weight mass, size and terminal velocity:
        real_t delta = IWCS / (IWCS + IWCL);
        real_t amass = delta * m_as + (real_t(1) - delta) * m_al;
        return std::max(real_t(1e-18), amass) * si::kilograms;
      }

      // mean terminal velocity of ice A particle
      // eq. A.15b in Grabowski 1999
      template <typename real_t>
      quantity<si::velocity, real_t> velocity_iceA(
        const quantity<si::dimensionless, real_t> &ria,
        const quantity<si::mass_density, real_t> &rhod_0
      ) {
        using namespace common;
        real_t IWC = std::max(real_t(1.e-9) * si::dimensionless(), rhod_0 * ria / si::kilograms * si::cubic_meters);
        real_t IWCS = std::max(real_t(1.e-9),
                               std::min(std::min(real_t(1.e-3), IWC),
                                        real_t(2.52e-4) * std::pow((IWC / real_t(1.e-3)), real_t(0.837))));
        real_t IWCL = std::max(real_t(1.e-9), IWC - IWCS);
        //velocity of small ice A particle:
        real_t v_as = real_t(0.1);
        //velocity of large ice A particle:
        real_t v_al = real_t(0.9) + real_t(0.1) * std::log10(real_t(1e3) * IWCL);
        //average velocity of ice A particle:
        real_t delta = IWCS / (IWCS + IWCL);
        return (delta * v_as + (real_t(1) - delta) * v_al) * std::pow(
          real_t(1.) / rhod_0 * si::kilograms / si::cubic_meters, real_t(0.5)) * si::meters / si::seconds;
      }


      // graupel density used for ice B (from Grabowski 1999)
      libcloudphxx_const(si::mass_density, rho_ib, 400, si::kilograms / si::cubic_metres)

      // slope of Marshall Palmer distribution for ice B
      // eq. A.4 in Grabowski 1999
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::length>::type, real_t> lambda_ice_b(
        const quantity<si::dimensionless, real_t> &rib,
        const quantity<si::mass_density, real_t> &rhod_0
      ) {
        using namespace common;
        real_t N_0b = real_t(4e6);
        // adding numerical epsilon to ensure there's no division by 0
        return std::pow(
          pi<real_t>() * rho_ib<real_t>() * N_0b / (rhod_0 * rib + std::numeric_limits<real_t>::epsilon() *
            si::kilograms / si::cubic_meters), real_t(0.25)) / si::meters;
      }

      // mean mass of ice B particle
      // eq. A.5 in Grabowski 1999
      template <typename real_t>
      quantity<si::mass, real_t> mass_b(
        const quantity<si::dimensionless, real_t> &rib,
        const quantity<si::mass_density, real_t> &rhod_0
      ) {
        using namespace common;
        real_t bmass = pi<real_t>() * rho_ib<real_t>() / (real_t(6) * std::pow(lambda_ice_b(rib, rhod_0) * si::meters,
                                                                               real_t(3))) * si::cubic_meters /
          si::kilograms;
        return std::max(real_t(1e-18), bmass) * si::kilograms;
      }

      // mean terminal velocity of ice B particle
      // eq. A.6 in Grabowski 1999
      template <typename real_t>
      quantity<si::velocity, real_t> velocity_iceB(
        const quantity<si::dimensionless, real_t> &rib,
        const quantity<si::mass_density, real_t> &rhod_0
      ) {
        return real_t(31.2) *
          std::pow(lambda_ice_b(rib, rhod_0) * si::meters, real_t(-0.37)) * std::pow(
            rhod_0 / si::kilograms * si::cubic_meters,
            real_t(-0.5)) *
          si::meters / si::seconds;
      }

      //coefficients for ice growth by deposition/riming
      // table 2 from Koenig 1972
      template <typename real_t>
      real_t coeff_alpha(const quantity<si::temperature, real_t> &T) {
        real_t alpha_list[32] = {
          real_t(0.), real_t(0.7939e-7), real_t(0.7841e-6), real_t(0.3369e-5), real_t(0.4336e-5), real_t(0.5285e-5),
          real_t(0.3728e-5), real_t(0.1852e-5), real_t(0.2991e-6), real_t(0.4248e-6), real_t(0.7434e-6),
          real_t(0.1812e-5),
          real_t(0.4394e-5), real_t(0.9145e-5), real_t(0.1725e-4), real_t(0.3348e-4), real_t(0.1725e-4),
          real_t(0.9175e-5),
          real_t(0.4412e-5), real_t(0.2252e-5), real_t(0.9115e-6), real_t(0.4876e-6), real_t(0.3473e-6),
          real_t(0.4758e-6),
          real_t(0.6306e-6), real_t(0.8573e-6), real_t(0.7868e-6), real_t(0.7192e-6), real_t(0.6515e-6),
          real_t(0.5956e-6),
          real_t(0.533e-6), real_t(0.4834e-6)
        };
        //temperature in Celsius:
        real_t Tc = T / si::kelvins - real_t(273.16);
        //limiting the temerature to (-31, 0):
        real_t ttcoe = std::min(real_t(0.), std::max(real_t(-31), Tc));
        int index1 = static_cast<int>(std::trunc(-ttcoe)); //index between 0 and 31
        int index2 = index1 + 1; // index between 1 and 32
        real_t del = -ttcoe - real_t(index1);
        return (real_t(1.) - del) * alpha_list[index1] + del * alpha_list[index2]; //interpolation
      }

      //coefficients for ice growth by deposition/riming
      // table 2 from Koenig 1972
      template <typename real_t>
      real_t coeff_beta(const quantity<si::temperature, real_t> &T) {
        real_t beta_list[32] = {
          real_t(0.), real_t(0.4006), real_t(0.4831), real_t(0.5320), real_t(0.5307), real_t(0.5319), real_t(0.5249),
          real_t(0.4888),
          real_t(0.3894), real_t(0.4047), real_t(0.4318), real_t(0.4771), real_t(0.5183), real_t(0.5463),
          real_t(0.5651), real_t(0.5813),
          real_t(0.5655), real_t(0.5478), real_t(0.5203), real_t(0.4906), real_t(0.4447), real_t(0.4126),
          real_t(0.3960), real_t(0.4149),
          real_t(0.4320), real_t(0.4506), real_t(0.4483), real_t(0.4460), real_t(0.4433), real_t(0.4413),
          real_t(0.4382), real_t(0.4361)
        };
        //temperature in Celsius:
        real_t Tc = T / si::kelvins - real_t(273.16);
        //limiting the temerature to (-31, 0):
        real_t ttcoe = std::min(real_t(0.), std::max(real_t(-31), Tc));
        int index1 = static_cast<int>(std::trunc(-ttcoe)); //index between 0 and 31
        int index2 = index1 + 1; // index between 1 and 32
        real_t del = -ttcoe - real_t(index1);
        return (real_t(1.) - del) * beta_list[index1] + del * beta_list[index2]; //interpolation
      }

      //ice A homogeneous nucleation 1 (rv -> ria)
      // eq. A.21a in Grabowski 1999
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> hom_A_nucleation_1(
        const quantity<si::dimensionless, real_t> &rv,
        const quantity<si::dimensionless, real_t> &rvs,
        const quantity<si::dimensionless, real_t> &rvsi,
        const quantity<si::temperature, real_t> &T,
        const quantity<si::time, real_t> &dt
      ) {
        const quantity<si::time, real_t> taunuc = dt; // timescale for nucleation
        real_t beta = (T > real_t(213.16) * si::kelvins)
                        ? real_t(0.1) + real_t(0.9) * (T - real_t(213.16) * si::kelvins) / (real_t(20) * si::kelvins)
                        : real_t(0.1);
        quantity<si::dimensionless, real_t> rv_adj = beta * rvs + (real_t(1) - beta) * rvsi;
        //homogeneous nucleation rate (only if T < -40 C)
        return (T < real_t(233.16) * si::kelvins)
                 ? (real_t(1) - std::exp(-dt / taunuc)) * std::max(real_t(0) * si::dimensionless(), rv - rv_adj) /
                 si::seconds
                 : real_t(0) / si::seconds;
      }

      //ice A homogeneous nucleation 2 (rc -> ria)
      // eq. A.21b in Grabowski 1999
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> hom_A_nucleation_2(
        const quantity<si::dimensionless, real_t> &rc,
        const quantity<si::temperature, real_t> &T,
        const quantity<si::time, real_t> &dt
      ) {
        const quantity<si::time, real_t> taunuc = dt; // time scale for nucleation
        //homogeneous nucleation rate (only if T < -40 C)
        return (T < real_t(233.16) * si::kelvins)
                 ? (real_t(1) - std::exp(-dt / taunuc)) * rc / si::seconds
                 : real_t(0) / si::seconds;
      }

      //ice A heterogeneous nucleation (rc-> ria)
      // eq. A.19 in Grabowski 1999
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> het_A_nucleation(
        const quantity<si::dimensionless, real_t> &ria,
        const quantity<si::dimensionless, real_t> &rc,
        const quantity<si::temperature, real_t> &T,
        const quantity<si::mass_density, real_t> &rhod_0,
        const quantity<si::time, real_t> &dt
      ) {
        if (rc == real_t(0) || T > real_t(273.16) * si::kelvins)
        {
          return real_t(0) / si::seconds;
        }
        using namespace common;
        const quantity<si::time, real_t> taunuc = dt; // timescale for nucleation
        //mean ice A particle mass
        quantity<si::mass, real_t> m_a = mass_a(ria, T, rhod_0);
        //concentration of ice nuclei:
        quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> N_in = std::min(real_t(1e5),
          real_t(1e-2) * std::exp(real_t(0.6) * (real_t(273.16) - T / si::kelvins))) / si::cubic_meters;
        //nucleation rate:
        real_t t_term = real_t(1) - std::exp(-dt / taunuc);
        return t_term * std::min(rc, std::max(real_t(0) * si::dimensionless(), (N_in * m_a / rhod_0) - ria)) /
          si::seconds;
      }

      // ice B heterogeneous nucleation 1 (rr -> rib)
      // eq. A.23 in Grabowski 1999
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> het_B_nucleation_1(
        const quantity<si::dimensionless, real_t> &rr,
        const quantity<si::dimensionless, real_t> &ria,
        const quantity<si::temperature, real_t> &T,
        const quantity<si::mass_density, real_t> &rhod_0
      ) {
        using namespace common;
        if (ria == real_t(0) || rr == real_t(0) || T > real_t(273.16) * si::kelvins)
        // ice B is formed only if ice A and rain are present
        {
          return real_t(0) / si::seconds;
        }
        quantity<divide_typeof_helper<si::dimensionless, si::length>::type, real_t> lambda_r = lambda_rain(rr, rhod_0);
        //mean raindrop mass (eq. A2 in Grabowski 1999)
        quantity<si::mass, real_t> m_r = pi<real_t>() * moist_air::rho_w<real_t>() / (real_t(6) * std::pow(
          lambda_r * si::meters, real_t(3))) * si::cubic_meters;
        //mean raindrop velocity (eq. A3 in Grabowski 1999)
        real_t v_r = real_t(251) * std::pow(lambda_r * si::meters * rhod_0 / si::kilograms * si::cubic_meters,
                                            real_t(-0.5));
        //mean raindrop radius (eq. A2 in Grabowski 1999)
        quantity<si::length, real_t> R_r = real_t(0.5) / lambda_r;
        //mean ice A particle mass
        quantity<si::mass, real_t> m_a = mass_a(ria, T, rhod_0);
        //mean ice A velocity
        real_t v_a = velocity_iceA(ria, rhod_0) * si::seconds / si::meters;
        //rate of collisions between raindrops and ice A
        real_t N_0r = real_t(1e7);
        quantity<divide_typeof_helper<si::frequency, si::mass>::type, real_t> N_ra = N_0r / si::cubic_meters /
          si::meters / lambda_r * std::abs(v_r - v_a) * si::meters / si::seconds * pi<real_t>() * R_r * R_r * ria / m_a;
        //nucleation rate:
        return N_ra * m_r;
      }

      //ice B heterogeneous nucleation 2 (ria -> rib)
      // eq. A.23 in Grabowski 1999
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> het_B_nucleation_2(
        const quantity<si::dimensionless, real_t> &rr,
        const quantity<si::dimensionless, real_t> &ria,
        const quantity<si::temperature, real_t> &T,
        const quantity<si::mass_density, real_t> &rhod_0
      ) {
        using namespace common;
        if (ria == real_t(0) || rr == real_t(0) || T > real_t(273.16) * si::kelvins)
        // ice B is formed only if ice A and rain are present
        {
          return real_t(0) / si::seconds;
        }
        quantity<divide_typeof_helper<si::dimensionless, si::length>::type, real_t> lambda_r = lambda_rain(rr, rhod_0);
        //mean raindrop velocity (eq. A3 in Grabowski 1999)
        real_t v_r = real_t(251) * std::pow(lambda_r * si::meters * rhod_0 / si::kilograms * si::cubic_meters,
                                            real_t(-0.5));
        //mean raindrop radius (eq. A2 in Grabowski 1999)
        quantity<si::length, real_t> R_r = real_t(0.5) / lambda_r;
        //mean ice A particle mass
        quantity<si::mass, real_t> m_a = mass_a(ria, T, rhod_0);
        //mean ice A velocity
        real_t v_a = velocity_iceA(ria, rhod_0) * si::seconds / si::meters;
        //rate of collisions between raindrops and ice A
        real_t N_0r = real_t(1e7);
        quantity<divide_typeof_helper<si::frequency, si::mass>::type, real_t> N_ra = N_0r / si::cubic_meters /
          si::meters / lambda_r * abs(v_r - v_a) * si::meters / si::seconds * pi<real_t>() * R_r * R_r * ria / m_a;
        //nucleation rate:
        return N_ra * m_a;
      }

      //melting of ice A (ria -> rr)
      // eq. A.26 in Grabowski 1999
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> melting_A(
        const quantity<si::dimensionless, real_t> &ria,
        const quantity<si::temperature, real_t> &T,
        const quantity<si::mass_density, real_t> &rhod_0,
        const quantity<si::time, real_t> &dt
      ) {
        using namespace common;
        if (ria == real_t(0) || T < real_t(273.16) * si::kelvin)
        // calculating melting rate only if ice A is present and T > 0 C
        {
          return real_t(0) / si::seconds;
        }
        //mass of average ice A particle:
        quantity<si::mass, real_t> m_a = mass_a(ria, T, rhod_0);
        //diameter of average ice A particle:
        quantity<si::length, real_t> D_a = std::pow(m_a / real_t(0.025) / si::kilograms, real_t(0.5)) * si::meters;
        //fall velocity of avg. ice A particle:
        quantity<si::velocity, real_t> v_a = velocity_iceA(ria, rhod_0);
        //Reynolds number:
        real_t Re = D_a * v_a * rhod_0 / vterm::visc(T);
        //ventilation factor:
        real_t F_a = std::max(real_t(1), real_t(0.78) + real_t(0.27) * std::pow(Re, real_t(0.5)));
        //melting rate:
        real_t dma_dt = real_t(9.e-7) * D_a / real_t(2) / si::meters * F_a * std::max(
          real_t(0) * si::dimensionless(), T / si::kelvins - real_t(273.16) * si::dimensionless());
        return std::min(ria / dt * si::seconds, dma_dt * ria / m_a * si::kilograms) / si::seconds;
      }

      //melting of ice B (rib -> rr)
      // eq. A.26 in Grabowski 1999
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> melting_B(
        const quantity<si::dimensionless, real_t> &rib,
        const quantity<si::temperature, real_t> &T,
        const quantity<si::mass_density, real_t> &rhod_0,
        const quantity<si::time, real_t> &dt
      ) {
        using namespace common;
        if (rib == real_t(0) || T < real_t(273.16) * si::kelvin)
        // calculating melting rate only if ice B is present and T > 0 C
        {
          return real_t(0) / si::seconds;
        }
        quantity<divide_typeof_helper<si::dimensionless, si::length>::type, real_t> lambda_b =
          lambda_ice_b(rib, rhod_0);
        //mass of average ice B particle:
        quantity<si::mass, real_t> m_b = mass_b(rib, rhod_0);
        //diameter of average ice B particle (eq. A5 in Grabowski 1999):
        quantity<si::length, real_t> D_b = real_t(1) / lambda_b;
        //fall velocity of avg. ice B particle:
        quantity<si::velocity, real_t> v_b = velocity_iceB(rib, rhod_0);
        //Reynolds number:
        real_t Re = D_b * v_b * rhod_0 / vterm::visc(T);
        //ventilation factor:
        real_t F_b = std::max(real_t(1), real_t(0.78) + real_t(0.27) * std::pow(Re, real_t(0.5)));
        //melting rate
        real_t dmb_dt = real_t(9.e-7) * D_b / real_t(2) / si::meters * F_b * std::max(
          real_t(0) * si::dimensionless(), T / si::kelvins - real_t(273.16) * si::dimensionless());
        return std::min(rib / dt * si::seconds, dmb_dt * rib / m_b * si::kilograms) / si::seconds;
      }


      // growth of ice A by deposition (rv->ria)
      // eq. A.24a in Grabowski 1999, eq. 27 in Koenig 1976
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> deposition_A(
        const quantity<si::dimensionless, real_t> &ria,
        const quantity<si::dimensionless, real_t> &rv,
        const quantity<si::dimensionless, real_t> &rvs,
        const quantity<si::dimensionless, real_t> &rvsi,
        const quantity<si::temperature, real_t> &T,
        const quantity<si::mass_density, real_t> &rhod_0
      ) {
        if (ria == real_t(0) || T > real_t(273.16) * si::kelvins) // growth only if ice A is present
        {
          return real_t(0) / si::seconds;
        }
        //mass of average ice A particle:
        quantity<si::mass, real_t> m_a = mass_a(ria, T, rhod_0);
        real_t alpha = coeff_alpha(T);
        real_t beta = coeff_beta(T);
        // growth rate for a single particle:
        quantity<divide_typeof_helper<si::mass, si::time>::type, real_t> dm_dt_AE = real_t(1e-3) * (rv - rvsi) / (rvs -
            rvsi + std::numeric_limits<real_t>::epsilon()) * alpha * // adding numerical epsilon to avoid division by 0
          std::pow(m_a * real_t(1e3) / si::kilograms, beta) * si::kilograms / si::seconds;
        //1e-3 comes from conversion g/sec into kg/sec
        return ria / m_a * dm_dt_AE;
      }


      //growth of ice A by riming (rc->ria)
      // eq. A.24b in Grabowski 1999, eq. 27-34 in Koenig 1976
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> riming_A(
        const quantity<si::dimensionless, real_t> &ria,
        const quantity<si::dimensionless, real_t> &rc,
        const quantity<si::dimensionless, real_t> &rv,
        const quantity<si::dimensionless, real_t> &rvs,
        const quantity<si::dimensionless, real_t> &rvsi,
        const quantity<si::temperature, real_t> &T,
        const quantity<si::mass_density, real_t> &rhod_0
      ) {
        if (ria == real_t(0) || rc == real_t(0) || T > real_t(273.16) * si::kelvins)
        // growth only if ice A and cloud droplets present
        {
          return real_t(0) / si::seconds;
        }
        //mass of average ice A particle:
        quantity<si::mass, real_t> m_a = mass_a(ria, T, rhod_0);
        real_t alpha = coeff_alpha(T);
        real_t beta = coeff_beta(T);
        // regime AE
        quantity<divide_typeof_helper<si::mass, si::time>::type, real_t> dm_dt_AE = real_t(1e-3) * (rv - rvsi) / (rvs -
            rvsi + std::numeric_limits<real_t>::epsilon()) * alpha * // adding numerical epsilon to avoid division by 0
          std::pow(m_a * real_t(1e3) / si::kilograms, beta) * si::kilograms / si::seconds;
        //1e-3 comes from conversion g/sec into kg/sec
        // regime BC
        real_t tan_theta = real_t(1.) + real_t(0.1) * std::log(
          rhod_0 * rc * real_t(1e3) / si::kilograms * si::cubic_meters);
        real_t gamma = alpha * std::pow(real_t(5e-8), beta);
        quantity<divide_typeof_helper<si::mass, si::time>::type, real_t> dm_dt_BC = real_t(1e-3) * gamma * std::pow(
          m_a / real_t(5e-11) / si::kilograms, tan_theta) * si::kilograms / si::seconds;
        //regime CD
        real_t dzeta = gamma * std::pow(real_t(2e3), tan_theta);
        real_t xi = std::log(rc * rhod_0 / si::kilograms * si::cubic_meters * real_t(1e9) / dzeta) / std::log(
          real_t(1e4));
        quantity<divide_typeof_helper<si::mass, si::time>::type, real_t> dm_dt_CD = real_t(1e-3) * dzeta * std::pow(
          m_a * real_t(1e7) / si::kilograms, xi) * si::kilograms / si::seconds;
        //total growth
        quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> rima = real_t(0) / si::seconds;
        if (m_a > real_t(5e-11) * si::kilograms && m_a <= real_t(1e-7) * si::kilograms)
        {
          rima += std::max(real_t(0) * si::kilograms / si::seconds, dm_dt_BC - dm_dt_AE) * ria / m_a;
        }
        if (m_a > real_t(1e-7) * si::kilograms)
        {
          rima += std::max(real_t(0) * si::kilograms / si::seconds, dm_dt_CD - dm_dt_AE) * ria / m_a;
        }
        return rima;
      }

      // growth of ice B by deposition (rv->rib)
      // eq. A.24c in Grabowski 1999, eq. 27 in Koenig 1976
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> deposition_B(
        const quantity<si::dimensionless, real_t> &rib,
        const quantity<si::dimensionless, real_t> &rv,
        const quantity<si::dimensionless, real_t> &rvs,
        const quantity<si::dimensionless, real_t> &rvsi,
        const quantity<si::temperature, real_t> &T,
        const quantity<si::mass_density, real_t> &rhod_0
      ) {
        if (rib == real_t(0) || T > real_t(273.16) * si::kelvins) // calculating growth rate only if ice B is present
        {
          return real_t(0) / si::seconds;
        }
        //mass of average ice B particle:
        quantity<si::mass, real_t> m_b = mass_b(rib, rhod_0);
        real_t alpha = coeff_alpha(T);
        real_t beta = coeff_beta(T);
        // growth rate for a single particle:
        quantity<divide_typeof_helper<si::mass, si::time>::type, real_t> dm_dt_AE = real_t(1e-3) * (rv - rvsi) / (rvs -
            rvsi + std::numeric_limits<real_t>::epsilon()) * alpha * // adding numerical epsilon to avoid division by 0
          std::pow(m_b * real_t(1e3) / si::kilograms, beta) * si::kilograms / si::seconds;
        //1e-3 comes from conversion g/sec into kg/sec
        return rib / m_b * dm_dt_AE;
      }

      //growth of ice B by riming  (rc, rr ->rib)
      // eq. A.24d in Grabowski 1999, eq. 27-34 in Koenig 1976
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> riming_B(
        const quantity<si::dimensionless, real_t> &rib,
        const quantity<si::dimensionless, real_t> &rc,
        const quantity<si::dimensionless, real_t> &rv,
        const quantity<si::dimensionless, real_t> &rvs,
        const quantity<si::dimensionless, real_t> &rvsi,
        const quantity<si::temperature, real_t> &T,
        const quantity<si::mass_density, real_t> &rhod_0
      ) {
        if (rib == real_t(0) || T > real_t(273.16) * si::kelvins) // calculating growth rate only if ice B is present
        {
          return real_t(0) / si::seconds;
        }
        //mass of average ice B particle:
        quantity<si::mass, real_t> m_b = mass_b(rib, rhod_0);
        real_t alpha = coeff_alpha(T);
        real_t beta = coeff_beta(T);
        // regime AE
        quantity<divide_typeof_helper<si::mass, si::time>::type, real_t> dm_dt_AE = real_t(1e-3) * (rv - rvsi) / (rvs -
            rvsi + std::numeric_limits<real_t>::epsilon()) * alpha * // adding numerical epsilon to avoid division by 0
          std::pow(m_b * real_t(1e3) / si::kilograms, beta) * si::kilograms / si::seconds;
        //1e-3 comes from conversion g/sec into kg/sec
        // regime BC
        real_t tan_theta = real_t(1.) + real_t(0.1) * std::log(
          rhod_0 * rc * real_t(1e3) / si::kilograms * si::cubic_meters);
        real_t gamma = alpha * std::pow(real_t(5e-8), beta);
        quantity<divide_typeof_helper<si::mass, si::time>::type, real_t> dm_dt_BC = real_t(1e-3) * gamma * std::pow(
          m_b / real_t(5e-11) / si::kilograms, tan_theta) * si::kilograms / si::seconds;
        //regime CD
        real_t dzeta = gamma * std::pow(real_t(2e3), tan_theta);
        real_t xi = std::log(rc * rhod_0 / si::kilograms * si::cubic_meters * real_t(1e9) / dzeta) / std::log(
          real_t(1e4));
        quantity<divide_typeof_helper<si::mass, si::time>::type, real_t> dm_dt_CD = real_t(1e-3) * dzeta * std::pow(
          m_b * real_t(1e7) / si::kilograms, xi) * si::kilograms / si::seconds;
        //total growth
        quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> rimb = real_t(0) / si::seconds;
        if (m_b > real_t(5e-11) * si::kilograms && m_b <= real_t(1e-7) * si::kilograms)
        {
          rimb += std::max(real_t(0) * si::kilograms / si::seconds, dm_dt_BC - dm_dt_AE) * rib / m_b;
        }
        if (m_b > real_t(1e-7) * si::kilograms)
        {
          rimb += std::max(real_t(0) * si::kilograms / si::seconds, dm_dt_CD - dm_dt_AE) * rib / m_b;
        }
        return rimb;
      }

      //growth of ice B by riming (only rc->rib)
      // eq. A.24d in Grabowski 1999, eq. 27-34 in Koenig 1976
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> riming_B_1(
        const quantity<si::dimensionless, real_t> &rib,
        const quantity<si::dimensionless, real_t> &rc,
        const quantity<si::dimensionless, real_t> &rr,
        const quantity<si::dimensionless, real_t> &rv,
        const quantity<si::dimensionless, real_t> &rvs,
        const quantity<si::dimensionless, real_t> &rvsi,
        const quantity<si::temperature, real_t> &T,
        const quantity<si::mass_density, real_t> &rhod_0
      ) {
        real_t coeff_rc = rc / (rc + rr + real_t(1e-10));
        return coeff_rc * riming_B(rib, rc, rv, rvs, rvsi, T, rhod_0);
      }

      //growth of ice B by riming (only rr->rib)
      // eq. A.24d in Grabowski 1999, eq. 27-34 in Koenig 1976
      template <typename real_t>
      quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> riming_B_2(
        const quantity<si::dimensionless, real_t> &rib,
        const quantity<si::dimensionless, real_t> &rc,
        const quantity<si::dimensionless, real_t> &rr,
        const quantity<si::dimensionless, real_t> &rv,
        const quantity<si::dimensionless, real_t> &rvs,
        const quantity<si::dimensionless, real_t> &rvsi,
        const quantity<si::temperature, real_t> &T,
        const quantity<si::mass_density, real_t> &rhod_0
      ) {
        real_t coeff_rc = rc / (rc + rr + real_t(1e-10));
        return (real_t(1) - coeff_rc) * riming_B(rib, rc, rv, rvs, rvsi, T, rhod_0);
      }
    };
  };
};
