/** @file
  * @copyright University of Warsaw
  * @brief Autoconversion and collection righ-hand side terms using Kessler formulae
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <algorithm>
#include <libcloudph++/common/detail/zip.hpp>
#include <libcloudph++/blk_2m/common_formulae.hpp>
#include <libcloudph++/blk_2m/activation_formulae.hpp>
#include <libcloudph++/blk_2m/cond_evap_formulae.hpp>

namespace libcloudphxx
{
  namespace blk_2m
  {
    template <typename real_t, class container_t>
    void forcings_elementwise(
      const opts_t<real_t> &opt,
      container_t drhod_th_cont,
      container_t drhod_rv_cont,
      container_t drhod_rc_cont,
      container_t drhod_nc_cont,
      const container_t rhod_cont,   
      const container_t rhod_th_cont,
      const container_t rhod_rv_cont,
      const container_t rhod_rc_cont,
      const container_t rhod_nc_cont
    )   
    {
      for (auto tup : zip(drhod_th_cont, drhod_rv_cont, drhod_rc_cont, drhod_nc_cont, 
                rhod_cont, rhod_th_cont,  rhod_rv_cont,  rhod_rc_cont,  rhod_nc_cont))
      {
        real_t
          &drhod_th = boost::get<0>(tup),
          &drhod_rv = boost::get<1>(tup),
          &drhod_rc = boost::get<2>(tup),
          &drhod_nc = boost::get<3>(tup);
        const quantity<si::mass_density, real_t> &rhod = boost::get<4>(tup) * si::kilograms / si::cubic_metres;
        const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t>  &rhod_th  = boost::get<5>(tup) * si::kilograms /si::cubic_metres * si::kelvin;
        const quantity<si::mass_density, real_t> &rhod_rv = boost::get<6>(tup) * si::kilograms / si::cubic_metres;
        const quantity<si::mass_density, real_t> &rhod_rc = boost::get<7>(tup) * si::kilograms / si::cubic_metres;
        const quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t>  &rhod_nc = boost::get<8>(tup) / si::cubic_metres;

        assert(opt.dt != 0);
        using namespace formulae;

        //helper temperature and pressure
        quantity<si::temperature, real_t> T = common::theta_dry::T<real_t>(rhod_th, rhod);
        quantity<si::pressure, real_t> p    = common::theta_dry::p<real_t>(rhod, rhod_rv/rhod, T);

        // if opt activation
        if (opt.acti)
        {
          if(rhod_rv/rhod > common::const_cp::r_vs<real_t>(T, p))
          {
            real_t tmp = activation_rate<real_t>(p, T, rhod, rhod_rv, rhod_nc, 
              opt.mean_rd, opt.sdev_rd, opt.N_tot, opt.dt * si::seconds, opt.chem_b) * si::seconds * si::cubic_metres; 

	    drhod_nc += tmp;  //TODO maybe make common ofr all the forcings?
            drhod_rv -= tmp * (ccnmass<real_t>() / si::kilograms);
            drhod_rc += tmp * (ccnmass<real_t>() / si::kilograms);
            drhod_th += tmp * (ccnmass<real_t>() / si::kilograms) 
                        * (common::theta_dry::d_rhodtheta_d_rv<real_t>(T, rhod_th) * si::cubic_metres / si::kilograms / si::kelvin); 

          }
        }

        // if opt condensation/evaporation
        if (opt.cond)
        {
          //TODO is it possible
          if(rhod_rc * si::cubic_metres / si::kilograms > 0 && rhod_nc * si::cubic_metres > 0)
          {
            real_t tmp = cond_evap_rate<real_t>(T, p, r_drop(rhod_rc, rhod_nc), rhod_rv/rhod, rhod_nc) * si::seconds
                         * (rhod * si::cubic_metres / si::kilograms);

            if(rhod_rc * si::cubic_metres / si::kilograms + tmp < 0) //so that we don't evaporate more than we have
            {
              real_t tmp = - rhod_rc * si::cubic_metres / si::kilograms;
            }

            drhod_rc += tmp;
            drhod_rv -= tmp;
            drhod_th -= tmp * (common::theta_dry::d_rhodtheta_d_rv<real_t>(T, rhod_th) * si::cubic_metres / si::kilograms / si::kelvin);

          }
        }
      }
    }
  }    
};
