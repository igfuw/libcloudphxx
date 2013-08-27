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
#include <libcloudph++/blk_2m/autoconversion_formulae.hpp>
#include <libcloudph++/blk_2m/accretion_formulae.hpp>
#include <libcloudph++/blk_2m/collision_sink_formulae.hpp>

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
      container_t drhod_rr_cont,
      container_t drhod_nr_cont,
      const container_t rhod_cont,   
      const container_t rhod_th_cont,
      const container_t rhod_rv_cont,
      const container_t rhod_rc_cont,
      const container_t rhod_nc_cont,
      const container_t rhod_rr_cont,
      const container_t rhod_nr_cont
    )   
    {  
      // sanity checks
      assert(opt.dt != 0);
      assert(min(rhod_rv_cont) > 0);
      assert(min(rhod_th_cont) > 0);
      assert(min(rhod_rc_cont) >= 0);
      assert(min(rhod_rr_cont) >= 0);
      assert(min(rhod_nc_cont) >= 0);
      assert(min(rhod_nr_cont) >= 0);

      assert(min(drhod_nc_cont) == 0);
      assert(min(drhod_nr_cont) == 0);
      assert(min(drhod_rc_cont) == 0);
      assert(min(drhod_rr_cont) == 0);
      assert(max(drhod_nc_cont) == 0);
      assert(max(drhod_nr_cont) == 0);
      assert(max(drhod_rc_cont) == 0);
      assert(max(drhod_rr_cont) == 0);

      using namespace formulae;
      using namespace common::moist_air;
      using namespace common::theta_dry;

      //timestep 
      quantity<si::time, real_t > dt = opt.dt * si::seconds;

      // if something is too small e-179 it becomes negative
      // so instead of rhod_nl == 0 we have rhod_nl < eps
      // also -1e-30 + 1e-30 is not equal to zero
      quantity<si::dimensionless, real_t>                                         eps_d = 1e-20;
      quantity<si::mass_density, real_t>                                          eps_r = 1e-20 * si::kilograms / si::cubic_metres;
      quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> eps_n = 1e-20 / si::cubic_metres;

                      //TODO: 
                      //unfortunately can't zip through more than 10 arguments 
                      //so instead one loop over all forcings, there will be a few 
      for (auto tup : zip(drhod_th_cont, drhod_rv_cont, drhod_rc_cont, drhod_nc_cont,
                rhod_cont, rhod_th_cont,  rhod_rv_cont,  rhod_rc_cont,  rhod_nc_cont))
      {
        real_t
          &drhod_th = boost::get<0>(tup),
          &drhod_rv = boost::get<1>(tup),
          &drhod_rc = boost::get<2>(tup),
          &drhod_nc = boost::get<3>(tup);
        const quantity<si::mass_density, real_t> &rhod    = boost::get<4>(tup) * si::kilograms / si::cubic_metres;
        const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t>  &rhod_th  = boost::get<5>(tup) * si::kilograms /si::cubic_metres * si::kelvin;
        const quantity<si::mass_density, real_t> &rhod_rv = boost::get<6>(tup) * si::kilograms / si::cubic_metres;
        const quantity<si::mass_density, real_t> &rhod_rc = boost::get<7>(tup) * si::kilograms / si::cubic_metres;
        const quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t>  &rhod_nc = boost::get<8>(tup) / si::cubic_metres;

        //helper temperature and pressure
        quantity<si::temperature, real_t> T = common::theta_dry::T<real_t>(rhod_th, rhod);
        quantity<si::pressure, real_t>    p = common::theta_dry::p<real_t>(rhod, rhod_rv/rhod, T);

        // activation (see Morrison & Grabowski 2007)
        if(opt.acti)
        { //TODO what if we have some oter source terms (that happen somewhere before here), like diffusion?
          assert(drhod_rc == 0 && "activation is first");
          assert(drhod_nc == 0 && "activation is first");
          assert(drhod_th == 0 && "activation is first");

          if(rhod_rv/rhod > common::const_cp::r_vs<real_t>(T, p))
          {
            quantity<divide_typeof_helper<si::frequency, si::volume>::type, real_t> tmp = 
              activation_rate<real_t>(p, T, rhod, rhod_rv, rhod_nc, opt.mean_rd, opt.sdev_rd, opt.N_stp, dt, opt.chem_b); 

	    drhod_nc += tmp * si::cubic_metres * si::seconds;  

            drhod_rv -= tmp * ccnmass<real_t>() / si::kilograms * si::cubic_metres * si::seconds;
            drhod_rc += tmp * ccnmass<real_t>() / si::kilograms * si::cubic_metres * si::seconds;

            //TODO maybe some common part for all the forcings (for example drhod_th)?
            drhod_th -= tmp * ccnmass<real_t>() * d_rhodtheta_d_rv<real_t>(T, rhod_th) / rhod
                         / si::kilograms / si::kelvins * si::cubic_metres * si::seconds; 
          }

          assert(drhod_nc >= 0 && "activation can only increase cloud droplet concentration");
          assert(drhod_rc >= 0 && "activation can only increase cloud water");
          assert(drhod_th >= 0 && "activation can only increase theta");
         }

        // condensation/evaporation of cloud water (see Morrison & Grabowski 2007)
        if(opt.cond)
        {                                       
          if (rhod_rc > eps_r && rhod_nc > eps_n)
          {                                              //  ^^   TODO is it possible?
            quantity<divide_typeof_helper<si::mass_density, si::time>::type, real_t> tmp = 
              cond_evap_rate<real_t>(T, p, r_drop(rhod_rc, rhod_nc), rhod_rv/rhod, rhod_nc) * rhod;

            if (rhod_rc + ((drhod_rc * si::kilograms / si::cubic_metres / si::seconds + tmp) * dt)  < 0 * si::kilograms / si::cubic_metres)             
            {   //so that we don't evaporate more cloud water than there is
              tmp      =- (rhod_rc + (dt * drhod_rc * si::kilograms / si::cubic_metres / si::seconds)) / dt;  //evaporate all rhod_rc
              drhod_nc =- rhod_nc / dt * si::cubic_metres * si::seconds; //and all rhod_nc
            }

            drhod_rc += tmp / si::kilograms * si::cubic_metres * si::seconds;
            drhod_rv -= tmp / si::kilograms * si::cubic_metres * si::seconds;

            drhod_th -= tmp  * d_rhodtheta_d_rv<real_t>(T, rhod_th) / rhod
                         / si::kilograms / si::kelvins * si::cubic_metres * si::seconds; 
          }

          assert(rhod_rc * si::cubic_metres / si::kilograms + drhod_rc * opt.dt >= 0 && "condensation/evaporation can't make rhod_rc < 0");
          assert(rhod_rv * si::cubic_metres / si::kilograms + drhod_rv * opt.dt >= 0 && "condensation/evaporation can't make rhod_rv < 0");
          assert(rhod_th * si::cubic_metres / si::kilograms / si::kelvin + drhod_th * opt.dt >= 0 && "condensation/evaporation can't make rhod_th < 0");
        }
      }

      for (auto tup : zip(drhod_rc_cont, drhod_nc_cont, drhod_rr_cont, drhod_nr_cont,
                rhod_cont, rhod_rc_cont,  rhod_nc_cont,  rhod_rr_cont))
      {
        real_t
          &drhod_rc = boost::get<0>(tup),
          &drhod_nc = boost::get<1>(tup),
          &drhod_rr = boost::get<2>(tup),
          &drhod_nr = boost::get<3>(tup);
        const quantity<si::mass_density, real_t> &rhod    = boost::get<4>(tup) * si::kilograms / si::cubic_metres;
        const quantity<si::mass_density, real_t> &rhod_rc = boost::get<5>(tup) * si::kilograms / si::cubic_metres;
        const quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t>  &rhod_nc = boost::get<6>(tup) / si::cubic_metres;
        const quantity<si::mass_density, real_t> &rhod_rr = boost::get<7>(tup) * si::kilograms / si::cubic_metres;
 
        if (rhod_rc * si::cubic_metres / si::kilograms + drhod_rc * opt.dt > 0)
        {
          //autoconversion rate (as in Khairoutdinov and Kogan 2000, but see Wood 2005 table 1)
          if(opt.acnv)
          {                                  
           if(rhod_rc > eps_r && rhod_nc > eps_n)
            {                     
              quantity<divide_typeof_helper<si::mass_density, si::time>::type, real_t> tmp = autoconv_rate(rhod, rhod_rc, rhod_nc);

              //so that autoconversion doesn't take more rhod_rc than there is
              tmp = std::min(tmp, (rhod_rc + dt * drhod_rc * si::kilograms / si::cubic_metres / si::seconds) / dt);
              assert(tmp * si::seconds * si::cubic_metres / si::kilograms >= 0 && "autoconv rate has to be >= 0");

              drhod_rc -= tmp / si::kilograms * si::cubic_metres * si::seconds;
              drhod_rr += tmp / si::kilograms * si::cubic_metres * si::seconds;

              //sink of N for cloud droplets is combined with the sink due to accretion
              //source of N for drizzle assumes that all the drops have the same radius
              drhod_nr += tmp / (4./3 * pi<real_t>() * pow<3>(drizzle_radius<real_t>()) * rho_w<real_t>())
                            * si::cubic_metres * si::seconds;
            }

            assert(rhod_rc * si::cubic_metres / si::kilograms + drhod_rc * opt.dt >= 0 && "autoconversion can't make rhod_rc negative");
          }

          if (rhod_rc * si::cubic_metres / si::kilograms + drhod_rc * opt.dt > 0)
          {
            //accretion rate (as in Khairoutdinov and Kogan 2000, but see Wood 2005 table 1)
            if(opt.accr)
            {              
              if (rhod_rc > eps_r && rhod_nc > eps_n && rhod_rr > eps_r)  
              {                   
                quantity<divide_typeof_helper<si::mass_density, si::time>::type, real_t> tmp = accretion_rate(rhod, rhod_rc, rhod_rr);
                //so that accretion doesn't take more rhod_rc than there is
                tmp = std::min(tmp, (rhod_rc + dt * drhod_rc * si::kilograms / si::cubic_metres / si::seconds) / dt);
                assert(tmp * si::seconds * si::cubic_metres / si::kilograms >= 0 && "accretion rate has to be >= 0");
          
                drhod_rr += tmp / si::kilograms * si::cubic_metres * si::seconds;
                drhod_rc -= tmp / si::kilograms * si::cubic_metres * si::seconds;
                //the sink of N for cloud droplets is combined with sink due to autoconversion
                //accretion does not change N for drizzle 
              }

              assert(rhod_rc * si::cubic_metres / si::kilograms + drhod_rc * opt.dt >= 0 && "accretion can't make rhod_rc negative");
            }
          }

          //sink of cloud droplet concentration due to autoconversion and accretion (see Khairoutdinov and Kogan 2000 eq 35)
          //                                                                        (be careful cause "q" there actually means mixing ratio)
          //has to be just after autoconv. and accretion so that drhod_rr is a sum of only those two
          if(opt.acnv || opt.accr)
          {
            if (rhod_nc > eps_n && drhod_rr > eps_d)  
            {                           
              quantity<divide_typeof_helper<si::frequency, si::volume>::type, real_t> tmp =
                collision_sink_rate(drhod_rr * si::kilograms / si::cubic_metres / si::seconds, r_drop(rhod_rc, rhod_nc));

              assert(tmp >= 0 / si::cubic_metres / si::seconds && "tmp");
 
              //so that collisions don't take more rhod_nc than there is
              tmp = std::min(tmp, (rhod_nc + dt * drhod_nc / si::cubic_metres / si::seconds) / dt);
 
              drhod_nc -= tmp * si::cubic_metres * si::seconds;
            }
          
          assert(rhod_nc * si::cubic_metres + drhod_nc * opt.dt >= 0 && "collisions can't make rhod_nc negative");
          } 
        }
      }

      for (auto tup : zip(drhod_th_cont, drhod_rv_cont, drhod_rr_cont, drhod_nr_cont,
                rhod_cont, rhod_th_cont,  rhod_rv_cont,  rhod_rr_cont,  rhod_nr_cont))
      {
        real_t
          &drhod_th = boost::get<0>(tup),
          &drhod_rv = boost::get<1>(tup),
          &drhod_rr = boost::get<2>(tup),
          &drhod_nr = boost::get<3>(tup);
        const quantity<si::mass_density, real_t> &rhod    = boost::get<4>(tup) * si::kilograms / si::cubic_metres;
        const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t>  &rhod_th  = boost::get<5>(tup) * si::kilograms /si::cubic_metres * si::kelvin;
        const quantity<si::mass_density, real_t> &rhod_rv = boost::get<6>(tup) * si::kilograms / si::cubic_metres;
        const quantity<si::mass_density, real_t> &rhod_rr = boost::get<7>(tup) * si::kilograms / si::cubic_metres;
        const quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t>  &rhod_nr = boost::get<8>(tup) / si::cubic_metres;

        //helper temperature and pressure
        quantity<si::temperature, real_t> T = common::theta_dry::T<real_t>(rhod_th, rhod);
        quantity<si::pressure, real_t>    p = common::theta_dry::p<real_t>(rhod, rhod_rv/rhod, T);

        // evaporation of rain (see Morrison & Grabowski 2007)
        if(opt.cond)
        {
          if(rhod_rr > eps_r && rhod_nr > eps_n)
          { //cond/evap for rhod_rr

            assert(rhod_rr * si::cubic_metres / si::kilograms + drhod_rr * opt.dt >= 0 && "before rain cond-evap");
            assert(rhod_rv * si::cubic_metres / si::kilograms + drhod_rv * opt.dt >= 0 && "before rain cond-evap");
            assert(rhod_nr * si::cubic_metres + drhod_nr * opt.dt >= 0 && "before rain cond-evap");
            assert(rhod_th * si::cubic_metres / si::kilograms / si::kelvin + drhod_th * opt.dt >= 0 && "before rain cond-evap");

            quantity<divide_typeof_helper<si::mass_density, si::time>::type, real_t> tmp = 
              cond_evap_rate<real_t>(T, p, r_drop(rhod_rr, rhod_nr), rhod_rv/rhod, rhod_nr) * rhod;
              //TODO ventilation coefficents (not so important in drizzle but very needed in rain)

            tmp=std::min(tmp , real_t(0) * si::kilograms / si::cubic_metres / si::seconds);
 
            if(rhod_rr + (drhod_rr * si::kilograms / si::cubic_metres / si::seconds + tmp) * dt < 0 * si::kilograms / si::cubic_metres) 
            //so that we don't evaporate more than we have
            {
              tmp = - (rhod_rr + dt * drhod_rr * si::kilograms / si::cubic_metres / si::seconds) / dt; //evaporate all rhod_rr

              drhod_rv -= tmp / si::kilograms * si::cubic_metres * si::seconds;
              drhod_rr += tmp / si::kilograms * si::cubic_metres * si::seconds;

              drhod_nr  = -rhod_nr / dt * si::cubic_metres * si::seconds; //and all rhod_nr
       
              drhod_th += -tmp  * d_rhodtheta_d_rv<real_t>(T, rhod_th) / rhod
                          / si::kilograms / si::kelvins * si::cubic_metres * si::seconds; 
            }
            else
            {
              drhod_rv -= tmp / si::kilograms * si::cubic_metres * si::seconds;
              drhod_rr += tmp / si::kilograms * si::cubic_metres * si::seconds;
           
              drhod_th += -tmp  * d_rhodtheta_d_rv<real_t>(T, rhod_th) / rhod
                         / si::kilograms / si::kelvins * si::cubic_metres * si::seconds; 
              //during evaporation rhod_nr is reduced so that a constant mean drizzle/raindrop radius is mantained
              if (tmp < 0 * si::kilograms / si::cubic_metres / si::seconds) 
              {
                quantity<divide_typeof_helper<si::frequency, si::volume>::type, real_t> drhod_nr_tmp = tmp * rhod_nr / rhod_rr;

                if(rhod_nr + (drhod_nr / si::cubic_metres / si::seconds + drhod_nr_tmp) * dt > 0 / si::cubic_metres)
                {
                  drhod_nr += drhod_nr_tmp * si::cubic_metres * si::seconds;
                }
               // else do nothing
              }
            }
          }

          assert(rhod_rr * si::cubic_metres / si::kilograms + drhod_rr * opt.dt >= 0 && "rain condensation/evaporation can't make rhod_rr < 0");
          assert(rhod_rv * si::cubic_metres / si::kilograms + drhod_rv * opt.dt >= 0 && "rain condensation/evaporation can't make rhod_rv < 0");
          assert(rhod_nr * si::cubic_metres + drhod_nr * opt.dt >= 0 && "rain condensation/evaporation can't make rhod_nr < 0");
          assert(rhod_th * si::cubic_metres / si::kilograms / si::kelvin + drhod_th * opt.dt >= 0 && "rain condensation/evaporation can't make rhod_th < 0");
        }
      }
    }
  }    
};
