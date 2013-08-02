/** @file
  * @copyright University of Warsaw
  * @brief Autoconversion and collection righ-hand side terms using Kessler formulae
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <algorithm>
#include <libcloudph++/common/detail/zip.hpp>
#include <libcloudph++/blk_2m/activation_formulae.hpp>

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

//std::cerr<<"rhod_theta = "<<rhod_th_cont<<std::endl;
//std::cerr<<"rhod_rv = "<<rhod_rv_cont<<std::endl;
//std::cerr<<"rhod_rc = "<<rhod_rc_cont<<std::endl;
//std::cerr<<"rhod_nc = "<<rhod_nc_cont<<std::endl;

      // activation
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

//std::cerr<<"beta default = "<< beta_default<real_t>()<<std::endl;
//std::cerr<<"ccnmass = "<< ccnmass<real_t>()<<std::endl;
//std::cerr<<"s_0 = "<<s_0(T, opt.mean_rd, opt.chem_b)<<std::endl;
//std::cerr<<"s = "<<s(p, T, rhod, rhod_rv)<<std::endl;
//std::cerr<<"sdev_rd_s = "<<sdev_rd_s(opt.sdev_rd)<<std::endl;
//std::cerr<<"u = "<<(p, T, rhod, rhod_rv, opt.mean_rd, opt.sdev_rd, opt.chem_b)<<std::endl;
//std::cerr<<"N_c_p = "<<N_c_p(p, T, rhod, rhod_rv, opt.mean_rd, opt.sdev_rd, opt.N_tot, opt.chem_b)<<std::endl;
//std::cerr<<"opt.dt = "<<opt.dt<<std::endl;
//std::cerr<<"activation rate = "<<activation_rate<real_t>(p, T, rhod, rhod_rv, rhod_nc, opt.mean_rd, opt.sdev_rd, opt.N_tot, opt.dt * si::seconds, opt.chem_b) * si::seconds * si::cubic_metres<<std::endl; 

//std::cerr<<"drhod_nc = "<< tmp << std::endl;
	    drhod_nc += tmp;
            drhod_rv -= tmp * (ccnmass<real_t>() / si::kilograms);
            drhod_rc += tmp * (ccnmass<real_t>() / si::kilograms);
            drhod_th += tmp * (ccnmass<real_t>() / si::kilograms) 
                        * (common::theta_dry::d_rhodtheta_d_rv<real_t>(T, rhod_th) * si::cubic_metres / si::kilograms / si::kelvin); 

          }
        }
      }
    }    
  }
};
