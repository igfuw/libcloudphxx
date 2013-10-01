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
//<listing>
    template <typename real_t, class cont_t>
    void rhs_cellwise(
      const opts_t<real_t> &opts,
      cont_t &dot_rho_e_cont,
      cont_t &dot_rho_v_cont,
      cont_t &dot_rho_c_cont,
      cont_t &dot_n_c_cont,
      cont_t &dot_rho_r_cont,
      cont_t &dot_n_r_cont,
      const cont_t &rho_d_cont,   
      const cont_t &rho_e_cont,
      const cont_t &rho_v_cont,
      const cont_t &rho_c_cont,
      const cont_t &n_c_cont,
      const cont_t &rho_r_cont,
      const cont_t &n_r_cont,
      const real_t &dt
    )   
//</listing>
    {  
      // sanity checks
      assert(min(rho_v_cont) > 0);
      assert(min(rho_e_cont) > 0);
      assert(min(rho_c_cont) >= 0);
      assert(min(rho_r_cont) >= 0);
      assert(min(n_c_cont) >= 0);
      assert(min(n_r_cont) >= 0);

      assert(min(dot_n_c_cont) == 0);
      //assert(min(dot_n_r_cont) == 0);
      assert(min(dot_rho_c_cont) == 0);
      //assert(min(dot_rho_r_cont) == 0);
      assert(max(dot_n_c_cont) == 0);
      assert(max(dot_n_r_cont) == 0);
      assert(max(dot_rho_c_cont) == 0);
      assert(max(dot_rho_r_cont) == 0);

      using namespace formulae;
      using namespace common::moist_air;
      using namespace common::theta_dry;

      // if something is too small e-179 it becomes negative
      // so instead of n_l == 0 we have n_l < eps
      // also -1e-30 + 1e-30 is not equal to zero
      quantity<si::dimensionless, real_t>                                         eps_d = real_t(1e-20);
      quantity<si::mass_density, real_t>                                          eps_r = real_t(1e-20) * si::kilograms / si::cubic_metres;
      quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> eps_n = real_t(1e-20) / si::cubic_metres;

                      //TODO: 
                      //unfortunately can't zip through more than 10 arguments 
                      //so instead one loop over all forcings, there will be a few 
      for (auto tup : zip(
        dot_rho_e_cont, 
        dot_rho_v_cont, 
        dot_rho_c_cont, 
        dot_n_c_cont,
        rho_d_cont, 
        rho_e_cont,  
        rho_v_cont,  
        rho_c_cont,  
        n_c_cont
      ))
      {
        real_t
          &dot_rho_e = boost::get<0>(tup),
          &dot_rho_v = boost::get<1>(tup),
          &dot_rho_c = boost::get<2>(tup),
          &dot_n_c = boost::get<3>(tup);
        const quantity<si::mass_density, real_t> &rho_d    = boost::get<4>(tup) * si::kilograms / si::cubic_metres;
        const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t>  &rho_e  = boost::get<5>(tup) * si::kilograms /si::cubic_metres * si::kelvin;
        const quantity<si::mass_density, real_t> &rho_v = boost::get<6>(tup) * si::kilograms / si::cubic_metres;
        const quantity<si::mass_density, real_t> &rho_c = boost::get<7>(tup) * si::kilograms / si::cubic_metres;
        const quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t>  &n_c = boost::get<8>(tup) / si::cubic_metres;

        //helper temperature and pressure
        quantity<si::temperature, real_t> T = common::theta_dry::T<real_t>(rho_e, rho_d);
        quantity<si::pressure, real_t>    p = common::theta_dry::p<real_t>(rho_d, rho_v / rho_d, T);

        // activation (see Morrison & Grabowski 2007)
        if(opts.acti)
        { //TODO what if we have some other source terms (that happen somewhere before here), like diffusion?
          assert(dot_rho_c == 0 && "activation is first");
          assert(dot_n_c == 0 && "activation is first");
          assert(dot_rho_e == 0 && "activation is first");

          if(rho_v / rho_d > common::const_cp::r_vs<real_t>(T, p))
          {
            quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t> N_ccn = 0;

            // looping over lognormal modes
            for (const auto &mode : opts.dry_distro)
            { 
              N_ccn += N_c_p<real_t>(
                p, T, rho_d, rho_v, 
                mode.mean_rd * si::metres, 
                mode.sdev_rd, 
                mode.N_stp / si::cubic_metres, 
                mode.chem_b,
                opts.RH_max
              ); 
            }

            // 
            quantity<divide_typeof_helper<si::frequency, si::volume>::type, real_t> tmp = 
              activation_rate<real_t>(N_ccn, n_c, dt * si::seconds);

	    dot_n_c += tmp * si::cubic_metres * si::seconds;  
            dot_rho_v -= tmp * ccnmass<real_t>() / si::kilograms * si::cubic_metres * si::seconds;
            dot_rho_c += tmp * ccnmass<real_t>() / si::kilograms * si::cubic_metres * si::seconds;

            //TODO maybe some common part for all the forcings (for example dot_rho_e)?
            dot_rho_e -= tmp * ccnmass<real_t>() * d_rhodtheta_d_rv<real_t>(T, rho_e) / rho_d
                         / si::kilograms / si::kelvins * si::cubic_metres * si::seconds; 
          }

          assert(dot_n_c >= 0 && "activation can only increase cloud droplet concentration");
          assert(dot_rho_c >= 0 && "activation can only increase cloud water");
          assert(dot_rho_e >= 0 && "activation can only increase theta");
         }

        // condensation/evaporation of cloud water (see Morrison & Grabowski 2007)
        if(opts.cond)
        {                                       
          if (rho_c > eps_r && n_c > eps_n)
          {             //  ^^   TODO is it possible?
            quantity<divide_typeof_helper<si::mass_density, si::time>::type, real_t> tmp = 
              cond_evap_rate<real_t>(T, p, rho_v / rho_d, tau_relax_c(T, p, r_drop_c(rho_c, n_c), n_c)) * rho_d;

std::cerr<<"promien kropelek "<<r_drop_c(rho_c, n_c)<<std::endl;

            assert(r_drop_c(rho_c, n_c) >= 0 * si::metres  && "mean droplet radius cannot be < 0");

            if (rho_c + ((dot_rho_c * si::kilograms / si::cubic_metres / si::seconds + tmp) * (dt * si::seconds))  < 0 * si::kilograms / si::cubic_metres)             
            {   //so that we don't evaporate more cloud water than there is
              tmp      =- (rho_c + (dt * dot_rho_c * si::kilograms / si::cubic_metres)) / (dt * si::seconds);  //evaporate all rho_c
              dot_n_c =- n_c / dt * si::cubic_metres; //and all rho_c
            }

            dot_rho_c += tmp / si::kilograms * si::cubic_metres * si::seconds;
            dot_rho_v -= tmp / si::kilograms * si::cubic_metres * si::seconds;

            dot_rho_e -= tmp  * d_rhodtheta_d_rv<real_t>(T, rho_e) / rho_d
                         / si::kilograms / si::kelvins * si::cubic_metres * si::seconds; 
          }
          assert(rho_c * si::cubic_metres / si::kilograms + dot_rho_c * dt >= 0 && "condensation/evaporation can't make rho_c < 0");
          assert(rho_v * si::cubic_metres / si::kilograms + dot_rho_v * dt >= 0 && "condensation/evaporation can't make rho_v < 0");
          assert(rho_e * si::cubic_metres / si::kilograms / si::kelvin + dot_rho_e * dt >= 0 && "condensation/evaporation can't make rho_e < 0");
        }
      }

      for (auto tup : zip(
        dot_rho_c_cont, 
        dot_n_c_cont, 
        dot_rho_r_cont, 
        dot_n_r_cont,
        rho_d_cont,
        rho_c_cont,  
        n_c_cont, 
        rho_r_cont
      ))
      {
        real_t
          &dot_rho_c = boost::get<0>(tup),
          &dot_n_c = boost::get<1>(tup),
          &dot_rho_r = boost::get<2>(tup),
          &dot_n_r = boost::get<3>(tup);
        const quantity<si::mass_density, real_t> &rho_d    = boost::get<4>(tup) * si::kilograms / si::cubic_metres;
        const quantity<si::mass_density, real_t> &rho_c = boost::get<5>(tup) * si::kilograms / si::cubic_metres;
        const quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t>  &n_c = boost::get<6>(tup) / si::cubic_metres;
        const quantity<si::mass_density, real_t> &rho_r = boost::get<7>(tup) * si::kilograms / si::cubic_metres;
 
        if (rho_c * si::cubic_metres / si::kilograms + dot_rho_c * dt > 0)
        {
          //autoconversion rate (as in Khairoutdinov and Kogan 2000, but see Wood 2005 table 1)
          if(opts.acnv)
          {                                  
           if(rho_c > eps_r && n_c > eps_n)
            {   
                  
              quantity<divide_typeof_helper<si::mass_density, si::time>::type, real_t> tmp = autoconv_rate(rho_d, rho_c, n_c);

              //so that autoconversion doesn't take more rho_c than there is
              tmp = std::min(tmp, (rho_c + dt * dot_rho_c * si::kilograms / si::cubic_metres) / (dt * si::seconds));
              assert(tmp * si::seconds * si::cubic_metres / si::kilograms >= 0 && "autoconv rate has to be >= 0");

              dot_rho_c -= tmp / si::kilograms * si::cubic_metres * si::seconds;
              dot_rho_r += tmp / si::kilograms * si::cubic_metres * si::seconds;

              //sink of N for cloud droplets is combined with the sink due to accretion
              //source of N for drizzle assumes that all the drops have the same radius
              dot_n_r += tmp / (4./3 * pi<real_t>() * pow<3>(drizzle_radius<real_t>()) * rho_w<real_t>())
                            * si::cubic_metres * si::seconds;
            }

            assert(rho_c * si::cubic_metres / si::kilograms + dot_rho_c * dt >= 0 && "autoconversion can't make rho_c negative");
          }

          if (rho_c * si::cubic_metres / si::kilograms + dot_rho_c * dt > 0)
          {
            //accretion rate (as in Khairoutdinov and Kogan 2000, but see Wood 2005 table 1)
            if(opts.accr)
            {              
              if (rho_c > eps_r && n_c > eps_n && rho_r > eps_r)  
              {                   
                quantity<divide_typeof_helper<si::mass_density, si::time>::type, real_t> tmp = accretion_rate(rho_d, rho_c, rho_r);
                //so that accretion doesn't take more rho_c than there is
                tmp = std::min(tmp, (rho_c + dt * dot_rho_c * si::kilograms / si::cubic_metres) / (dt * si::seconds));
                assert(tmp * si::seconds * si::cubic_metres / si::kilograms >= 0 && "accretion rate has to be >= 0");
          
                dot_rho_r += tmp / si::kilograms * si::cubic_metres * si::seconds;
                dot_rho_c -= tmp / si::kilograms * si::cubic_metres * si::seconds;
                //the sink of N for cloud droplets is combined with sink due to autoconversion
                //accretion does not change N for drizzle 
              }

              assert(rho_c * si::cubic_metres / si::kilograms + dot_rho_c * dt >= 0 && "accretion can't make rho_c negative");
            }
          }

          //sink of cloud droplet concentration due to autoconversion and accretion (see Khairoutdinov and Kogan 2000 eq 35)
          //                                                                        (be careful cause "q" there actually means mixing ratio)
          //has to be just after autoconv. and accretion so that dot_rho_r is a sum of only those two
          if(opts.acnv || opts.accr)
          {
            if (n_c > eps_n && dot_rho_r > eps_d)  
            {                           
              quantity<divide_typeof_helper<si::frequency, si::volume>::type, real_t> tmp =
                collision_sink_rate(dot_rho_r * si::kilograms / si::cubic_metres / si::seconds, r_drop_c(rho_c, n_c));

              assert(r_drop_c(rho_c, n_c) >= 0 * si::metres  && "mean droplet radius cannot be < 0");
              assert(tmp >= 0 / si::cubic_metres / si::seconds && "tmp");
 
              //so that collisions don't take more n_c than there is
              tmp = std::min(tmp, (n_c + dt * dot_n_c / si::cubic_metres) / (dt * si::seconds));
 
              dot_n_c -= tmp * si::cubic_metres * si::seconds;
            }
          
          assert(n_c * si::cubic_metres + dot_n_c * dt >= 0 && "collisions can't make n_c negative");
          } 
        }
      }

      for (auto tup : zip(
        dot_rho_e_cont,
        dot_rho_v_cont, 
        dot_rho_r_cont, 
        dot_n_r_cont,
        rho_d_cont,
        rho_e_cont, 
        rho_v_cont, 
        rho_r_cont,
        n_r_cont
      )) {
        real_t
          &dot_rho_e = boost::get<0>(tup),
          &dot_rho_v = boost::get<1>(tup),
          &dot_rho_r = boost::get<2>(tup),
          &dot_n_r = boost::get<3>(tup);
        const quantity<si::mass_density, real_t> &rho_d    = boost::get<4>(tup) * si::kilograms / si::cubic_metres;
        const quantity<multiply_typeof_helper<si::mass_density, si::temperature>::type, real_t>  &rho_e  = boost::get<5>(tup) * si::kilograms /si::cubic_metres * si::kelvin;
        const quantity<si::mass_density, real_t> &rho_v = boost::get<6>(tup) * si::kilograms / si::cubic_metres;
        const quantity<si::mass_density, real_t> &rho_r = boost::get<7>(tup) * si::kilograms / si::cubic_metres;
        const quantity<divide_typeof_helper<si::dimensionless, si::volume>::type, real_t>  &n_r = boost::get<8>(tup) / si::cubic_metres;

        //helper temperature and pressure
        quantity<si::temperature, real_t> T = common::theta_dry::T<real_t>(rho_e, rho_d);
        quantity<si::pressure, real_t>    p = common::theta_dry::p<real_t>(rho_d, rho_v / rho_d, T);

        // evaporation of rain (see Morrison & Grabowski 2007)
        if(opts.cond)
        {
//std::cerr << "before evap: rho_r = " << rho_r * si::cubic_metres / si::kilograms <<std::endl;
          if(rho_r > eps_r && n_r > eps_n)
          { //cond/evap for rho_r

            assert(rho_r * si::cubic_metres / si::kilograms + dot_rho_r * dt >= 0 && "before rain cond-evap");
            assert(rho_v * si::cubic_metres / si::kilograms + dot_rho_v * dt >= 0 && "before rain cond-evap");
            assert(n_r * si::cubic_metres + dot_n_r * dt >= 0 && "before rain cond-evap");
            assert(rho_e * si::cubic_metres / si::kilograms / si::kelvin + dot_rho_e * dt >= 0 && "before rain cond-evap");

            quantity<divide_typeof_helper<si::mass_density, si::time>::type, real_t> tmp = 
              cond_evap_rate<real_t>(T, p, rho_v / rho_d, tau_relax_r(T, rho_d, rho_r, n_r)) * rho_d;

            assert(r_drop_r(rho_r, n_r) >= 0 * si::metres  && "mean drop radius cannot be < 0");

std::cerr<<"promień kropel "<<r_drop_r(rho_r, n_r)<<std::endl;
std::cerr<<"  "<<std::endl;

            tmp=std::min(tmp , real_t(0) * si::kilograms / si::cubic_metres / si::seconds);
 
            if(rho_r + (dot_rho_r * si::kilograms / si::cubic_metres / si::seconds + tmp) * (dt * si::seconds) < 0 * si::kilograms / si::cubic_metres) 
            //so that we don't evaporate more than we have
            {
              tmp = - (rho_r + dt * dot_rho_r * si::kilograms / si::cubic_metres) / (dt * si::seconds); //evaporate all rho_r

              dot_rho_v -= tmp / si::kilograms * si::cubic_metres * si::seconds;
              dot_rho_r += tmp / si::kilograms * si::cubic_metres * si::seconds;

              dot_n_r  = -n_r / dt * si::cubic_metres; //and all n_r
       
              dot_rho_e += -tmp  * d_rhodtheta_d_rv<real_t>(T, rho_e) / rho_d
                          / si::kilograms / si::kelvins * si::cubic_metres * si::seconds; 
            }
            else
            {
              dot_rho_v -= tmp / si::kilograms * si::cubic_metres * si::seconds;
              dot_rho_r += tmp / si::kilograms * si::cubic_metres * si::seconds;
           
              dot_rho_e += -tmp  * d_rhodtheta_d_rv<real_t>(T, rho_e) / rho_d
                         / si::kilograms / si::kelvins * si::cubic_metres * si::seconds; 
              //during evaporation n_r is reduced so that a constant mean drizzle/raindrop radius is mantained
              if (tmp < 0 * si::kilograms / si::cubic_metres / si::seconds) 
              {
                quantity<divide_typeof_helper<si::frequency, si::volume>::type, real_t> dot_n_r_tmp = tmp * n_r / rho_r;

                if(n_r + (dot_n_r / si::cubic_metres / si::seconds + dot_n_r_tmp) * (dt * si::seconds) > 0 / si::cubic_metres)
                {
                  dot_n_r += dot_n_r_tmp * si::cubic_metres * si::seconds;
                }
               // else do nothing
              }
            }
          }

//std::cerr << "after evap rho_r = " << rho_r << std::endl;
//std::cerr << "dot_rho_r * dt = " << dot_rho_r * dt <<std::endl;
//std::cerr << "rho_r + dot_rho_r * dt " << rho_r * si::cubic_metres / si::kilograms + dot_rho_r * dt << std::endl; 

          assert(rho_r * si::cubic_metres / si::kilograms + dot_rho_r * dt >= 0 && "rain condensation/evaporation can't make rho_r < 0");
          assert(rho_v * si::cubic_metres / si::kilograms + dot_rho_v * dt >= 0 && "rain condensation/evaporation can't make rho_v < 0");
          assert(n_r * si::cubic_metres + dot_n_r * dt >= 0 && "rain condensation/evaporation can't make n_r < 0");
          assert(rho_e * si::cubic_metres / si::kilograms / si::kelvin + dot_rho_e * dt >= 0 && "rain condensation/evaporation can't make rho_e < 0");
        }
      }
    }
  };    
};
