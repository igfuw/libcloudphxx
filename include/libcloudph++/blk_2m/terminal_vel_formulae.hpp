/** @file
  * @copyright University of Warsaw
  * @brief double-moment bulk condensation/evaporation parameterisation formulae
  * @date August 2013
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once
#include <libcloudph++/common/moist_air.hpp>
#include <libcloudph++/blk_2m/common_formulae.hpp>

namespace libcloudphxx
{
  namespace blk_2m
  {
    namespace formulae
    { 
      //terminal fall velocity based on the data from Gunn and Kinzer (1949) and Beard (1976) 
      //modified by Simmel et al. (2002) -> see table 2 there

      libcloudphxx_const(si::length, d1, 134.43  * 1e-6, si::metres)
      libcloudphxx_const(si::length, d2, 1511.64 * 1e-6, si::metres)
      libcloudphxx_const(si::length, d3, 3477.84 * 1e-6, si::metres)

      template <typename real_t>
      quantity<si::dimensionless, real_t> alpha_fall(
        const quantity<si::mass_density, real_t> &rhod_rr,
        const quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> &rhod_nr
      ) {
         assert(r_drop_r(rhod_rr, rhod_nr) >= 0 * si::metres  && "mean drop radius cannot be < 0");
         
         if (real_t(2) * r_drop_r(rhod_rr, rhod_nr) == 0 * si::metres)   {return 0;} //TODO if no rain, terminal velocity = 0
         else if (real_t(2) * r_drop_r(rhod_rr, rhod_nr) < d1<real_t>()) {return 4.5795 * 1e5;} 
         else if (real_t(2) * r_drop_r(rhod_rr, rhod_nr) < d2<real_t>()) {return 4.962  * 1e3;} 
         else if (real_t(2) * r_drop_r(rhod_rr, rhod_nr) < d3<real_t>()) {return 1.732  * 1e3;} 
         else                                                            {return 9.17   * 1e2;} 
      }

      template <typename real_t>
      quantity<si::dimensionless, real_t> beta_fall(
        const quantity<si::mass_density, real_t> &rhod_rr,
        const quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> &rhod_nr
      ) {
         assert(r_drop_r(rhod_rr, rhod_nr) >= 0 * si::metres  && "mean drop radius cannot be < 0");

         if (real_t(2) * r_drop_r(rhod_rr, rhod_nr) < d1<real_t>())      {return 2./3;} 
         else if (real_t(2) * r_drop_r(rhod_rr, rhod_nr) < d2<real_t>()) {return 1./3;} 
         else if (real_t(2) * r_drop_r(rhod_rr, rhod_nr) < d3<real_t>()) {return 1./6;} 
         else                                                            {return 0;} 
      }

      template <typename real_t>
      quantity<si::velocity, real_t> v_term(
        const quantity<si::mass_density, real_t> &rhod_rr,
        const quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> &rhod_nr
      ) { 
            return alpha_fall(rhod_rr, rhod_nr)             //.... mass of drop in grams
                   * pow(rhod_rr/rhod_nr/si::kilograms * 1000, beta_fall(rhod_rr, rhod_nr)) 
                                        //^^^^^^^^^^^^ to make it dimensionless       
                   * real_t(1e-2) * si::metres/si::seconds; //velocity in metres/seconds
      }
    };
  };
};
