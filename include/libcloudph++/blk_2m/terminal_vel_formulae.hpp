/** @file
  * @copyright University of Warsaw
  * @brief double-moment bulk condensation/evaporation parameterisation formulae
  * @date August 2013
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once
#include <libcloudph++/common/moist_air.hpp>
#include <libcloudph++/common/earth.hpp>
#include <libcloudph++/blk_2m/common_formulae.hpp>

namespace libcloudphxx
{
  namespace blk_2m
  {
    namespace formulae
    {
      // terminal fall velocity based on the data from Gunn and Kinzer (1949) and Beard (1976)
      // modified by Simmel et al. (2002) -> see table 2 there

      // the actuall fall velocity for rain density or concentration is calculated as mass or number weighted mean:
      // v_m = \int N(D) m(D) vt(D) dD
      // v_n = \int N(D)      vt(D) dD

      libcloudphxx_const(si::length, d1,  134.43 * 1e-6, si::metres)
      libcloudphxx_const(si::length, d2, 1511.64 * 1e-6, si::metres)
      libcloudphxx_const(si::length, d3, 3477.84 * 1e-6, si::metres)

      template <typename real_t>
      inline quantity<si::dimensionless, real_t> alpha_fall(
        const quantity<si::length, real_t> &drop_r
      ) {
         assert(drop_r >= 0 * si::metres  && "mean drop radius cannot be < 0");

         if (real_t(2) * drop_r == 0 * si::metres)   {return 0;} //TODO if no rain, terminal velocity = 0
         else if (real_t(2) * drop_r < d1<real_t>()) {return 4.5795 * 1e5;}
         else if (real_t(2) * drop_r < d2<real_t>()) {return 4.962  * 1e3;}
         else if (real_t(2) * drop_r < d3<real_t>()) {return 1.732  * 1e3;}
         else                                        {return 9.17   * 1e2;}
      }

      template <typename real_t>
      inline quantity<si::dimensionless, real_t> beta_fall(
        const quantity<si::length, real_t> &drop_r
      ) {
         assert(drop_r >= 0 * si::metres  && "mean drop radius cannot be < 0");

         if (real_t(2) * drop_r < d1<real_t>())      {return 2./3;}
         else if (real_t(2) * drop_r < d2<real_t>()) {return 1./3;}
         else if (real_t(2) * drop_r < d3<real_t>()) {return 1./6;}
         else                                        {return 0;}
      }

      //helper intergrals for vm
        //int( D^5 * exp(-lbd * D) )
        template <typename real_t>
        inline quantity<si::dimensionless, real_t >mint_1(
          const quantity<divide_typeof_helper<si::dimensionless, si::length>::type, real_t> &lbd, //slope of assumed exponential size distribution
          const quantity<si::length, real_t> &D // diameter
        ) {
          auto tmp =  - std::pow(lbd * si::metres, -6) * exp(-lbd * D) *
            (pow(real_t(lbd * D), 5) + 5 * std::pow(real_t(lbd * D), 4) * 20 * std::pow(real_t(lbd * D), 3) + 60 * std::pow(real_t(lbd * D), 2) + 120 * real_t(lbd * D) + 120);

          assert(finite(tmp) && "mint_1 is finite failed");
          return tmp;
        }
        //int( D^4 * exp(-lbd * D) )
        template <typename real_t>
        inline quantity<si::dimensionless, real_t >mint_2(
          const quantity<divide_typeof_helper<si::dimensionless, si::length>::type, real_t> &lbd, //slope of assumed exponential size distribution
          const quantity<si::length, real_t> &D // diameter
        ) {
          auto tmp = - std::pow(real_t(lbd * si::metres), -5) * exp(-lbd * D) *
            (pow(real_t(lbd * D), 4) + 4 * std::pow(real_t(lbd * D), 3) * 12 * std::pow(real_t(lbd * D), 2) + 24 * real_t(lbd * D) + 24);

          assert(finite(tmp) && "mint_2 is finite failed");
          return tmp;
        }
        //int( D^(7/2) * exp(-lbd * D) )
        template <typename real_t>
        inline quantity<si::dimensionless, real_t >mint_3(
          const quantity<divide_typeof_helper<si::dimensionless, si::length>::type, real_t> &lbd, // slope of assumed exponential size distribution
          const quantity<si::length, real_t> &D // diameter
        ) {
          auto tmp = real_t(1./16) / std::pow(real_t(lbd * si::metres), real_t(9./2))
                   * (105 * std::sqrt(pi<real_t>()) * std::erf(sqrt(lbd * D))
                      - 2 * std::sqrt(lbd * D) * exp(-lbd * D) * (8 * std::pow(real_t(lbd * D), 3) + 28 * std::pow(real_t(lbd * D), 2) + 70 * real_t(lbd * D) + 105)
                   );

          assert(finite(tmp) && "mint_3 is finite failed");
          return tmp;
        }

        // int( exp(-lbd * D) )
        template <typename real_t>
        inline quantity<si::dimensionless, real_t >int_4(
          const quantity<divide_typeof_helper<si::dimensionless, si::length>::type, real_t> &lbd, //slope of assumed exponential size distribution
          const quantity<si::length, real_t> &D // diameter
        ) {
          auto tmp = - std::pow(real_t(lbd * si::metres), -1) * exp(-lbd * D);

          assert(finite(tmp) && "int_4 is finite failed");
          return tmp;
        }

       // helper intergrals for vn
       // int(D^(2) * exp(-lbd * D))
       template <typename real_t>
       inline quantity<si::dimensionless, real_t >nint_1(
         const quantity<divide_typeof_helper<si::dimensionless, si::length>::type, real_t> &lbd, // slope of assumed exponential size distribution
         const quantity<si::length, real_t> &D // diameter
       ) {
         auto tmp = std::pow(real_t(lbd * si::metres), -3) * exp(-lbd * D) * (- real_t(lbd * D) * (real_t(lbd * D) + 2) - 2);

         assert(finite(tmp) && "nint_1 is finite failed");
         return tmp;
       }

       // int(D * exp(-lbd * D))
       template <typename real_t>
       inline quantity<si::dimensionless, real_t >nint_2(
         const quantity<divide_typeof_helper<si::dimensionless, si::length>::type, real_t> &lbd, //slope of assumed exponential size distribution
         const quantity<si::length, real_t> &D // diameter
       ) {
         auto tmp =  - std::pow(real_t(lbd * si::metres), -2) * exp(-lbd * D) * (lbd * D + 1);

         assert(finite(tmp) && "nint_2 is finite failed");
         return tmp;
       }

       // int(D^(1/2) * exp(-lbd * D))
       template <typename real_t>
       inline quantity<si::dimensionless, real_t >nint_3(
         const quantity<divide_typeof_helper<si::dimensionless, si::length>::type, real_t> &lbd, //slope of assumed exponential size distribution
         const quantity<si::length, real_t> &D //diameter
       ) {
         auto tmp = std::sqrt(pi<real_t>()) * std::erf(sqrt(lbd * D)) / 2 / std::pow(real_t(lbd * si::metres), 3./2)
           - std::sqrt(D / si::metres) * exp(-lbd*D) / (lbd * si::metres);

         assert(finite(tmp) && "nint_3 is finite failed");
         return tmp;
       }

      template <typename real_t>
      inline quantity<si::velocity, real_t> v_term_m(
        const quantity<si::mass_density, real_t> &rhod,
        const quantity<si::dimensionless, real_t> &rr,
        const quantity<divide_typeof_helper<si::dimensionless, si::mass>::type, real_t> &nr
      ) {
        if (rr < rr_eps<real_t>() || nr < nr_eps<real_t>() / si::kilograms)
          return 0 * si::metres_per_second;

        quantity<divide_typeof_helper<si::dimensionless, si::length>::type, real_t> lbd = lambda_r(nr, rr);

        // eq. A4 in Morrison 2005
        auto tmp = rho_stp<real_t>() / rhod
          * lbd * si::metres * c_md<real_t>()
          * si::cubic_metres / si::kilograms // to make it dimensionless
          * real_t(1000) // mass of the drop in grams
          * (
            alpha_fall(d1<real_t>() / real_t(2))
            * std::pow(c_md<real_t>() * si::cubic_metres / si::kilograms * real_t(1000), beta_fall(d1<real_t>() / real_t(2)))
            * (mint_1(lbd, d1<real_t>()) - mint_1(lbd, real_t(0) * si::metres))
            +
            alpha_fall(d1<real_t>() + d2<real_t>() / real_t(2))
            * std::pow(c_md<real_t>() * si::cubic_metres / si::kilograms * real_t(1000), beta_fall(d1<real_t>() + d2<real_t>() / real_t(2)))
            * (mint_2(lbd, d2<real_t>()) - mint_2(lbd, d1<real_t>()))
            +
            alpha_fall(d2<real_t>() + d3<real_t>() / real_t(2))
            * std::pow(c_md<real_t>() * si::cubic_metres / si::kilograms * real_t(1000), beta_fall(d2<real_t>() + d3<real_t>() / real_t(2)))
            * (mint_3(lbd, d3<real_t>()) - mint_3(lbd, d2<real_t>()))
            +
            alpha_fall(real_t(2) * d3<real_t>())
            * (real_t(0) - int_4(lbd, d3<real_t>()))
          ) * real_t(1e-2) * si::metres/si::seconds;  // velocity in metres/seconds

//return std::max(real_t(0), real_t(tmp / si::metres_per_second)) * si::metres_per_second;
        assert(tmp * si::seconds / si::metres >= 0 && "negative terminal velocity!");
        assert(finite(tmp * si::seconds / si::metres) && "v_term_m terminal velocity is finite failed");
        return tmp;
      }

      template <typename real_t>
      inline quantity<si::velocity, real_t> v_term_n(
        const quantity<si::mass_density, real_t> &rhod,
        const quantity<si::dimensionless, real_t> &rr,
        const quantity<divide_typeof_helper<si::dimensionless, si::mass>::type, real_t> &nr
      ) {
        if (rr < rr_eps<real_t>() || nr < nr_eps<real_t>() / si::kilograms)
          return 0 * si::metres_per_second;

        quantity<divide_typeof_helper<si::dimensionless, si::length>::type, real_t> lbd = lambda_r(nr, rr);

        // eq A4 in Morrison 2005
        auto tmp = rho_stp<real_t>() / rhod
          * (
           alpha_fall(d1<real_t>() / real_t(2))
           * std::pow(c_md<real_t>() * si::cubic_metres / si::kilograms * real_t(1000), beta_fall(d1<real_t>() / real_t(2)))
           * (nint_1(lbd, d1<real_t>()) - nint_1(lbd, real_t(0) * si::metres))
           +
           alpha_fall(d1<real_t>() + d2<real_t>() / real_t(2))
           * std::pow(c_md<real_t>() * si::cubic_metres / si::kilograms * real_t(1000), beta_fall(d1<real_t>() + d2<real_t>() / real_t(2)))
           * (nint_2(lbd, d2<real_t>()) - nint_2(lbd, d1<real_t>()))
           +
           alpha_fall(d2<real_t>() + d3<real_t>() / real_t(2))
           * std::pow(c_md<real_t>() * si::cubic_metres / si::kilograms * real_t(1000), beta_fall(d2<real_t>() + d3<real_t>() / real_t(2)))
           * (nint_3(lbd, d3<real_t>()) - nint_3(lbd, d2<real_t>()))
           +
           alpha_fall(real_t(2) * d3<real_t>()) * (real_t(0) - int_4(lbd, d3<real_t>()))
          )
          * real_t(1e-2) * si::metres/si::seconds;  // velocity in metres/seconds

//return std::max(real_t(0), real_t(tmp / si::metres_per_second)) * si::metres_per_second;
        assert(tmp * si::seconds / si::metres >= 0 && "negative terminal velocity!");
        assert(finite(tmp * si::seconds / si::metres) && "v_term_n terminal velocity is finite failed");
        return tmp;
      }
    };
  };
};
