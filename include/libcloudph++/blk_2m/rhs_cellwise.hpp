/** @file
  * @copyright University of Warsaw
  * @brief Autoconversion and collection righ-hand side terms using Kessler formulae
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libcloudph++/blk_2m/extincl.hpp>

#include <libcloudph++/common/theta_dry.hpp>

namespace libcloudphxx
{
  namespace blk_2m
  {
//<listing>
    template <typename real_t, class cont_t>
    void rhs_cellwise(
      const opts_t<real_t> &opts,
      cont_t &dot_th_cont,
      cont_t &dot_rv_cont,
      cont_t &dot_rc_cont,
      cont_t &dot_nc_cont,
      cont_t &dot_rr_cont,
      cont_t &dot_nr_cont,
      const cont_t &rhod_cont,   
      const cont_t &th_cont,
      const cont_t &rv_cont,
      cont_t &rc_cont,
      cont_t &nc_cont,
      cont_t &rr_cont,
      cont_t &nr_cont,
      const real_t &dt
    )   
//</listing>
    { 
      real_t rc_eps = 1e-9;
      real_t nc_eps = 1e3;
      real_t rr_eps = 1e-12;
      real_t nr_eps = 1.;
/*
      for (auto tup : zip(
        rc_cont,
        nc_cont,
        rr_cont,
        nr_cont
      )) {
        real_t
          &rc = boost::get<0>(tup),
          &nc = boost::get<1>(tup),
          &rr = boost::get<2>(tup),
          &nr = boost::get<3>(tup);

        if (rc < rc_eps || nc < nc_eps)
        {
          rc = real_t(0);
          nc = real_t(0);
        }
        if (rr < rr_eps || nr < nr_eps)
        {
          rr = real_t(0);
          nr = real_t(0);
        }
      }
*/
      // sanity checks
      assert(min(rv_cont) >= 0);
      assert(min(th_cont) > 0);
      assert(min(rc_cont) >= 0);
      assert(min(rr_cont) >= 0);
      assert(min(nc_cont) >= 0);
      assert(min(nr_cont) >= 0); //<----------------------

      // TODO: rewrite so thet's not needed
      assert(min(dot_nc_cont) == 0);
      //assert(min(dot_nr_cont) == 0);
      assert(min(dot_rc_cont) == 0);
      //assert(min(dot_rr_cont) == 0);
      assert(max(dot_nc_cont) == 0);
      assert(max(dot_nr_cont) == 0);
      assert(max(dot_rc_cont) == 0);
      assert(max(dot_rr_cont) == 0);

      using namespace formulae;
      using namespace common::moist_air;
      using namespace common::theta_dry;

      //unfortunately can't zip through more than 10 arguments 
      //so instead one loop over all forcings, there will be a few 
      for (auto tup : zip(
        dot_th_cont, 
        dot_rv_cont, 
        dot_rc_cont, 
        dot_nc_cont,
        rhod_cont, 
        th_cont,  
        rv_cont,  
        rc_cont,  
        nc_cont
      ))
      {
        real_t
          &dot_th = boost::get<0>(tup),
          &dot_rv = boost::get<1>(tup),
          &dot_rc = boost::get<2>(tup),
          &dot_nc = boost::get<3>(tup);
        const quantity<si::mass_density,  real_t> &rhod  = boost::get<4>(tup) * si::kilograms / si::cubic_metres;
        const quantity<si::temperature,   real_t> &th    = boost::get<5>(tup) * si::kelvins;
        const quantity<si::dimensionless, real_t> &rv    = boost::get<6>(tup) * si::dimensionless();
        real_t &rc = boost::get<7>(tup);
        real_t &nc = boost::get<8>(tup);

        //helper temperature and pressure
        quantity<si::temperature, real_t> T = common::theta_dry::T<real_t>(th, rhod);
        quantity<si::pressure, real_t>    p = common::theta_dry::p<real_t>(rhod, rv, T);

        // activation (see Morrison & Grabowski 2007)
        if (opts.acti)
        { //TODO what if we have some other source terms (that happen somewhere before here), like diffusion?
          assert(dot_rc == 0 && "activation is first");
          assert(dot_nc == 0 && "activation is first");
          //assert(dot_th == 0 && "activation is first");

          if (rv > common::const_cp::r_vs<real_t>(T, p))
          {
            // summing by looping over lognormal modes
            quantity<divide_typeof_helper<si::dimensionless, si::mass>::type, real_t> n_ccn = 0;
            for (const auto &mode : opts.dry_distros)
            { 
              n_ccn += n_c_p<real_t>(
                p, T, rv, 
                mode.mean_rd * si::metres, 
                mode.sdev_rd, 
                mode.N_stp / si::cubic_metres, 
                mode.chem_b,
                opts.RH_max
              ); 
            }

            quantity<divide_typeof_helper<si::frequency, si::mass>::type, real_t> tmp = 
              activation_rate<real_t>(n_ccn, nc / si::kilograms, dt * si::seconds);

	    dot_nc += tmp * si::kilograms * si::seconds;  
            dot_rv -= tmp * ccnmass<real_t>() * si::seconds;
            dot_rc += tmp * ccnmass<real_t>() * si::seconds;

            //TODO maybe some common part for all the forcings (for example dot_th)?
            dot_th -= tmp * ccnmass<real_t>() * d_th_d_rv<real_t>(T, th) / si::kelvins * si::seconds; 
          }

          assert(dot_nc >= 0 && "activation can only increase cloud droplet concentration");
          assert(dot_rc >= 0 && "activation can only increase cloud water");
          //assert(dot_th >= 0 && "activation can only increase theta");
        }

        // condensation/evaporation of cloud water (see Morrison & Grabowski 2007)
        if (opts.cond)
        {                          
          if (rc > rc_eps && nc > nc_eps)
          {      //  ^^   TODO is it possible?
            quantity<divide_typeof_helper<si::dimensionless, si::time>::type, real_t> tmp = 
              cond_evap_rate<real_t>(
                T, p, rv, tau_relax_c(T, p, r_drop_c(rc, nc, rhod), rhod * nc / si::kilograms)
              );

            assert(r_drop_c(rc, nc, rhod) >= 0 * si::metres  && "mean droplet radius cannot be < 0");

            if (rc + ((dot_rc / si::seconds + tmp) * (dt * si::seconds))  < 0)
            {   // so that we don't evaporate more cloud water than there is
              tmp     = -(rc + (dt * dot_rc)) / (dt * si::seconds);  // evaporate all rc
              dot_nc  = -nc / dt;                                    // and all nc
            }
            dot_rc += tmp * si::seconds;

            if (rc + dot_rc * dt < 0)
            {  // (*)
                 // if rc is very small due to numerical reasons the above condition 
                 // may result in small negative values 
                 // (for rc = 1e-8 and dot_rc * dt = 1e-8 the new rc was -1e-25)
                 // in this case shamelessly put rc and nc to zero and do not calculate the impact on theta
std::cout << "SHAME!" << std::endl;
                tmp = 0;
                rc = 0;
                nc = 0;
                dot_rc = 0;
                dot_nc = 0;
            }
            dot_rv -= tmp * si::seconds;
            dot_th -= tmp  * d_th_d_rv<real_t>(T, th) / si::kelvins * si::seconds; 
          }

          assert(rc + dot_rc * dt >= 0 && "condensation/evaporation can't make rc < 0");
          assert(nc + dot_nc * dt >= 0 && "condensation/evaporation can't make nc < 0");
          assert(rv + dot_rv * dt >= 0 && "condensation/evaporation can't make rv < 0");
          assert(th / si::kelvin + dot_th * dt >= 0 && "condensation/evaporation can't make th < 0");
        }
      }

      for (auto tup : zip(
        dot_rc_cont, 
        dot_nc_cont, 
        dot_rr_cont, 
        dot_nr_cont,
        rhod_cont,
        rc_cont,  
        nc_cont, 
        rr_cont
      ))
      {
        real_t
          &dot_rc = boost::get<0>(tup),
          &dot_nc = boost::get<1>(tup),
          &dot_rr = boost::get<2>(tup),
          &dot_nr = boost::get<3>(tup);
        const quantity<si::mass_density, real_t>  &rhod = boost::get<4>(tup) * si::kilograms / si::cubic_metres;
        real_t &rc = boost::get<5>(tup);
        real_t &nc = boost::get<6>(tup);
        const quantity<si::dimensionless, real_t>   &rr = boost::get<7>(tup) * si::dimensionless();

        // autoconversion rate (as in Khairoutdinov and Kogan 2000, but see Wood 2005 table 1)
        if (opts.acnv)
        { 
          if (rc > rc_eps && nc > nc_eps)
          {  
            quantity<si::frequency, real_t> tmp = autoconv_rate(rc, nc, rhod, 
                                                                opts.acnv_A * si::dimensionless(), 
                                                                opts.acnv_b * si::dimensionless(), 
                                                                opts.acnv_c * si::dimensionless()
                                                               );

            // so that autoconversion doesn't take more rc than there is
            tmp = std::min(tmp, (rc + dt * dot_rc) / (dt * si::seconds));
            assert(tmp * si::seconds >= 0 && "autoconv rate has to be >= 0");

            dot_rc -= tmp * si::seconds;
            if (rc + dot_rc * dt < 0)
            { //see comment (*) in condensation 
              tmp = 0;
              rc = 0;
              dot_rc = 0;
            }
            dot_rr += tmp * si::seconds;

            // sink of N for cloud droplets is combined with the sink due to accretion
            // source of N for drizzle assumes that all the drops have the same radius
            dot_nr += tmp / (real_t(4)/3 * pi<real_t>() * rho_w<real_t>() * pow<3>(drizzle_radius<real_t>()))
              * si::kilograms * si::seconds; // to make it dimensionless
          }
          assert(rc + dot_rc * dt >= 0 && "autoconversion can't make rc negative");
        }

        // accretion rate (as in Khairoutdinov and Kogan 2000, but see Wood 2005 table 1)
        if (opts.accr)
        {              
          if (rc > rc_eps && nc > nc_eps && rr > rr_eps)  
          {                   
            quantity<si::frequency, real_t> tmp = accretion_rate(rc, rr);
            // so that accretion doesn't take more rc than there is
            tmp = std::min(tmp, (rc + dt * dot_rc) / (dt * si::seconds));
            assert(tmp * si::seconds >= 0 && "accretion rate has to be >= 0");
          
            dot_rc -= tmp * si::seconds;
            if (rc + dot_rc * dt < 0)
            { //see comment (*) in condensation 
              tmp = 0;
              rc = 0;
              dot_rc = 0;
            }
            dot_rr += tmp * si::seconds;

            // the sink of N for cloud droplets is combined with sink due to autoconversion
            // accretion does not change N for drizzle 
          }

          assert(rc + dot_rc * dt >= 0 && "accretion can't make rc negative");
        }

        // sink of n_c due to autoconversion and accretion (see Khairoutdinov and Kogan 2000 eq 35)
        //                                                 (be careful cause "q" there actually means mixing ratio, not water content)
        // has to be just after autoconv. and accretion so that dot_rr is a sum of only those two
        if (opts.acnv || opts.accr)
        {
          if (nc > nc_eps && dot_rr > 0)  
          {                           
            quantity<divide_typeof_helper<si::frequency, si::mass>::type, real_t> tmp =
              collision_sink_rate(dot_rr / si::seconds, r_drop_c(rc, nc, rhod));

            assert(r_drop_c(rc, nc, rhod) >= 0 * si::metres  && "mean droplet radius cannot be < 0");
            assert(tmp >= 0 / si::kilograms / si::seconds && "tmp");
 
            // so that collisions don't take more n_c than there is
            tmp = std::min(tmp, (nc / si::kilograms + dt * dot_nc / si::kilograms) / (dt * si::seconds));
            dot_nc -= tmp * si::kilograms * si::seconds;
            if (nc + dot_nc * dt < 0)
            { //see comment (*) in condensation
              nc = 0;
              dot_nc = 0;
            }
          }
          assert(nc + dot_nc * dt >= 0 && "collisions can't make n_c negative");
        } 
      }

      for (auto tup : zip(
        dot_th_cont,
        dot_rv_cont, 
        dot_rr_cont, 
        dot_nr_cont,
        rhod_cont,
        th_cont, 
        rv_cont, 
        rr_cont,
        nr_cont
      )) {
        real_t
          &dot_th = boost::get<0>(tup),
          &dot_rv = boost::get<1>(tup),
          &dot_rr = boost::get<2>(tup),
          &dot_nr = boost::get<3>(tup);
        const quantity<si::mass_density,  real_t> &rhod  = boost::get<4>(tup) * si::kilograms / si::cubic_metres;
        const quantity<si::temperature,   real_t> &th    = boost::get<5>(tup) * si::kelvins;
        const quantity<si::dimensionless, real_t> &rv    = boost::get<6>(tup) * si::dimensionless();
        real_t &rr    = boost::get<7>(tup);
        real_t &nr    = boost::get<8>(tup);

        quantity<si::dimensionless, real_t> rr_dim = rr * si::dimensionless();
        quantity<divide_typeof_helper<si::dimensionless, si::mass>::type, real_t> nr_dim = nr / si::kilograms;

        // helper temperature and pressure (TODO: it is repeated above!)
        quantity<si::temperature, real_t> T = common::theta_dry::T<real_t>(th, rhod);
        quantity<si::pressure, real_t>    p = common::theta_dry::p<real_t>(rhod, rv, T);

        // evaporation of rain (see Morrison & Grabowski 2007)
        if (opts.cond)
        {
          if (rr > rr_eps && nr_dim * si::kilograms > nr_eps)
          { // cond/evap for rr
            assert(rr_dim + dot_rr * dt >= 0 && "before rain cond-evap");
            assert(rv + dot_rv * dt >= 0 && "before rain cond-evap");
            assert(nr_dim * si::kilograms + dot_nr * dt >= 0 && "before rain cond-evap");
            assert(th / si::kelvin + dot_th * dt >= 0 && "before rain cond-evap");

            quantity<si::frequency, real_t> tmp = 
              cond_evap_rate<real_t>(T, p, rv, tau_relax_r(T, rhod, rr_dim, nr_dim));

            assert(r_drop_r(rr_dim, nr_dim) >= 0 * si::metres  && "mean drop radius cannot be < 0");

            tmp = std::min(tmp, real_t(0) / si::seconds);
            if (rr_dim + (dot_rr / si::seconds + tmp) * (dt * si::seconds) < 0) // so that we don't evaporate more than we have
            {
              tmp = - (rr_dim + dt * dot_rr) / (dt * si::seconds); // evaporate all rr
              dot_rv -= tmp * si::seconds;
              dot_th += -tmp  * d_th_d_rv<real_t>(T, th) / si::kelvins * si::seconds; 

              //dot_rr += tmp * si::seconds;
              //dot_nr  = -nr_dim / dt * si::kilograms; // and all n_r
              dot_rr = 0;
              dot_nr = 0;
              nr = 0;
              rr = 0;
            }
            else
            {
              dot_rv -= tmp * si::seconds;
              dot_rr += tmp * si::seconds;
           
              dot_th += -tmp * d_th_d_rv<real_t>(T, th) / si::kelvins * si::seconds; 

              // during evaporation n_r is reduced so that a constant mean drizzle/raindrop radius is mantained
              if (tmp < 0 / si::seconds) 
              {
                quantity<divide_typeof_helper<si::frequency, si::mass>::type, real_t> dot_nr_tmp = tmp * nr_dim / rr_dim;

                if (nr_dim + (dot_nr / si::kilograms / si::seconds + dot_nr_tmp) * (dt * si::seconds) > 0 / si::kilograms)
                {
                  dot_nr += dot_nr_tmp * si::kilograms * si::seconds;
                }
                // else do nothing
              }
            }
          }

          assert(rr_dim + dot_rr * dt >= 0 && "rain condensation/evaporation can't make rr < 0");
          assert(rv + dot_rv * dt >= 0 && "rain condensation/evaporation can't make rv < 0");
          assert(nr_dim * si::kilograms + dot_nr * dt >= 0 && "rain condensation/evaporation can't make n_r < 0");
          assert(th / si::kelvin + dot_th * dt >= 0 && "rain condensation/evaporation can't make re < 0");
        }
      }
/*
      //TODO - cleanup to be consistent with the eps from different microphysics source terms
      for (auto tup : zip(
        rc_cont,
        nc_cont,
        rr_cont,
        nr_cont
      )) {
        real_t
          &rc = boost::get<0>(tup),
          &nc = boost::get<1>(tup),
          &rr = boost::get<2>(tup),
          &nr = boost::get<3>(tup);

        if (rc < rc_eps || nc < nc_eps)
        {
          rc = real_t(0);
          nc = real_t(0);
        }
        if (rr < rr_eps || nr < nr_eps)
        {
          rr = real_t(0);
          nr = real_t(0);
        }
      }
*/
    }
  };    
};
