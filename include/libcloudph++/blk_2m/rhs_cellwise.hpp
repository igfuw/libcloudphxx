/** @file
  * @copyright University of Warsaw
  * @brief Autoconversion and collection righ-hand side terms using Kessler formulae
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libcloudph++/blk_2m/extincl.hpp>

#include <libcloudph++/common/theta_dry.hpp>
#include <libcloudph++/common/theta_std.hpp>

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
      const cont_t &th_cont,              // dry potential temperature (if const_p == false) or "standard" potential temperature (if const_p == true)
      const cont_t &rv_cont,
      const cont_t &rc_cont,
      const cont_t &nc_cont,
      const cont_t &rr_cont,
      const cont_t &nr_cont,
      const real_t &dt,
      const bool const_p = false,         // is pressure constant (e.g. anelastic approximation on UWLCM)?
      const cont_t &p_cont = cont_t()     // pressure, required if const_p == true
    )
//</listing>
    {
      // sanity checks
      assert(min(rv_cont) >= 0);
      assert(min(th_cont) > 0);
      assert(min(rc_cont) >= 0);
      assert(min(rr_cont) >= 0);
      assert(min(nc_cont) >= 0);
      assert(min(nr_cont) >= 0);
      assert(!const_p || p_cont.size() == th_cont.size());
      assert(!const_p || min(p_cont) > 0);

      using namespace formulae;
      using namespace common::moist_air;
      using namespace common::theta_dry;
      using namespace common::theta_std;

      for (auto tup : zip(
        dot_th_cont,
        dot_rv_cont,
        dot_rc_cont,
        dot_nc_cont,
        dot_rr_cont,
        dot_nr_cont,
        rhod_cont,
        th_cont,
        rv_cont,
        rc_cont,
        nc_cont,
        rr_cont,
        nr_cont,
        p_cont           // NOTE: for const_p == false, is zipping an empty p_cont safe?
      ))
      {
        real_t
          &dot_th = std::get<0>(tup),
          &dot_rv = std::get<1>(tup),
          &dot_rc = std::get<2>(tup),
          &dot_nc = std::get<3>(tup),
          &dot_rr = std::get<4>(tup),
          &dot_nr = std::get<5>(tup);
        const quantity<si::mass_density,  real_t> &rhod  = std::get<6>(tup) * si::kilograms / si::cubic_metres;
        const quantity<si::temperature,   real_t> &th    = std::get<7>(tup) * si::kelvins;
        const quantity<si::dimensionless, real_t> &rv    = std::get<8>(tup) * si::dimensionless();
        const real_t
          &rc = std::get<9>(tup),
          &nc = std::get<10>(tup),
          &rr = std::get<11>(tup),
          &nr = std::get<12>(tup);

        // helper dimensionless verions of real_t...
        const quantity<si::dimensionless, real_t> rr_dim = rr * si::dimensionless();
        // ... and a dimensional version of concentration
        const quantity<divide_typeof_helper<si::dimensionless, si::mass>::type, real_t> nr_dim = nr / si::kilograms;

        //helper temperature and pressure
        quantity<si::temperature, real_t> T;
        quantity<si::pressure, real_t>    p;

        if(!const_p)
        {
          T = common::theta_dry::T<real_t>(th, rhod);
          p = common::theta_dry::p<real_t>(rhod, rv, T);
        }
        else
        {
          p = std::get<13>(tup) * si::pascals;
          T = th * common::theta_std::exner(p);
        }

        // rhs only due to rhs_cellwise microphysics functions (needed for limiting)
        real_t local_dot_rc = 0,
               local_dot_rr = 0,
               local_dot_nc = 0,
               local_dot_nr = 0;

        // helper flags for when all cloud/rain evaporated and we skip collisions
        bool
          cloud_limiter = false,
          rain_limiter  = false;

        // activation (see Morrison & Grabowski 2007)
        if (opts.acti)
        {
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

	    local_dot_nc += tmp * si::kilograms * si::seconds;
            local_dot_rc += tmp * ccnmass<real_t>() * si::seconds;
          }

          assert(local_dot_nc >= 0 && "activation can only increase cloud droplet concentration");
          assert(local_dot_rc >= 0 && "activation can only increase cloud water");
        }

        if (opts.cond)
        {
          // condensation/evaporation of cloud water (see Morrison & Grabowski 2007)
          if (rc > rc_eps<real_t>() && nc > nc_eps<real_t>())
          {      //  ^^   TODO is it possible?
            quantity<si::frequency, real_t> tmp =
              cond_evap_rate<real_t>(
                T, p, rv, tau_relax_c(T, p, r_drop_c(rc, nc, rhod), rhod * nc / si::kilograms)
              );

            assert(r_drop_c(rc, nc, rhod) >= 0 * si::metres  && "mean droplet radius cannot be < 0");

            local_dot_rc += tmp * si::seconds;
          }

          // evaporation of rain (see Morrison & Grabowski 2007)
          if (rr > rr_eps<real_t>() && nr > nr_eps<real_t>())
          {
            quantity<si::frequency, real_t> tmp =
              std::min(
                cond_evap_rate<real_t>(T, p, rv, tau_relax_r(T, rhod, rr_dim, nr_dim)),
                real_t(0) / si::seconds  // only evaporation for rain
              );

            assert(r_drop_r(rr_dim, nr_dim) >= 0 * si::metres  && "mean drop radius cannot be < 0");

            local_dot_rr += tmp * si::seconds;
            // during evaporation n_r is reduced so that a constant mean drizzle/raindrop radius is mantained
            local_dot_nr += tmp * nr / rr * si::seconds;
          }
        }

        //limiters after acnv, cond/evap for rc and evap for rr
        local_dot_rc = std::max(local_dot_rc, - rc / dt);
        local_dot_rr = std::max(local_dot_rr, - rr / dt);
        local_dot_nr = std::max(local_dot_nr, - nr / dt);

        if(local_dot_rc == - rc / dt)
        {
          local_dot_nc = - nc / dt;
          cloud_limiter = true;
        }

        if(local_dot_rr == - rr / dt)
        {
          local_dot_nr = - nr / dt;
          rain_limiter  = true;
        }

        dot_rv -= (local_dot_rc + local_dot_rr);
        if(!const_p)
          dot_th -= (local_dot_rc + local_dot_rr) * d_th_d_rv<real_t>(T, th) / si::kelvins;
        else
          dot_th += common::const_cp::l_v(T) / (common::moist_air::c_pd<real_t>() * common::theta_std::exner(p)) * (local_dot_rc + local_dot_rr) / si::kelvins;
        dot_rc += local_dot_rc;
        dot_rr += local_dot_rr;
        dot_nc += local_dot_nc;
        dot_nr += local_dot_nr;

        // zero-out for reuse in collisions
        local_dot_rc = 0;
        local_dot_rr = 0;
        local_dot_nc = 0;
        local_dot_nr = 0;

        if (!cloud_limiter) // only do collisions if not all cloud water was evaporated
        {
          // autoconversion rate (as in Khairoutdinov and Kogan 2000, but see Wood 2005 table 1)
          if (opts.acnv)
          {
            if (rc > rc_eps<real_t>() && nc > nc_eps<real_t>())
            {
              quantity<si::frequency, real_t> tmp = autoconv_rate(rc, nc, rhod,
                                                                  opts.acnv_A * si::dimensionless(),
                                                                  opts.acnv_b * si::dimensionless(),
                                                                  opts.acnv_c * si::dimensionless()
                                                                 );

              // so that autoconversion doesn't take more rc than there is
              tmp = std::min(tmp, rc / (dt * si::seconds));

              assert(tmp * si::seconds >= 0 && "autoconv rate has to be >= 0");

              local_dot_rc -= tmp * si::seconds;
              local_dot_rr += tmp * si::seconds;

              // sink of N for cloud droplets is combined with the sink due to accretion
              // source of N for drizzle assumes that all the drops have the same radius
              local_dot_nr += tmp / (real_t(4)/3 * pi<real_t>() * rho_w<real_t>() *
                drizzle_radius<real_t>() * drizzle_radius<real_t>() * drizzle_radius<real_t>())
                * si::kilograms * si::seconds; // to make it dimensionless

              if(tmp == rc / (dt * si::seconds)) // all cloud water turned into rain
                cloud_limiter = true;
            }
          }

          // accretion rate (as in Khairoutdinov and Kogan 2000, but see Wood 2005 table 1)
          if (opts.accr && !cloud_limiter && !rain_limiter)
          {
            if (rc > rc_eps<real_t>() && nc > nc_eps<real_t>() && rr > rr_eps<real_t>())
            {
              quantity<si::frequency, real_t> tmp = accretion_rate(rc, rr_dim);

              assert(tmp * si::seconds >= 0 && "accretion rate has to be >= 0");

              local_dot_rc -= tmp * si::seconds;
              local_dot_rr += tmp * si::seconds;

              // limit dot_rc coming from acnv and accr
              local_dot_rc = std::max(local_dot_rc, - rc / dt);

              if(local_dot_rc == -rc / dt) // all cloud water turned into rain
                cloud_limiter = true;

              // the sink of N for cloud droplets is combined with sink due to autoconversion
              // accretion does not change N for drizzle
            }
          }

          // sink of n_c due to autoconversion and accretion
          //    (see Khairoutdinov and Kogan 2000 eq 35
          //     but be careful cause "q" there actually means mixing ratio, not water content)
          if (opts.acnv || opts.accr)
          {
            // if all cloud water was turned into rain, set nc = 0
            if(cloud_limiter)
              local_dot_nc = - nc / dt;
            // else calc sink of nc from local_dot_rr
            else if (nc > nc_eps<real_t>() && local_dot_rr > rr_eps<real_t>())
            {
              quantity<divide_typeof_helper<si::frequency, si::mass>::type, real_t> tmp =
                collision_sink_rate(local_dot_rr / si::seconds, r_drop_c(rc, nc, rhod));

              assert(r_drop_c(rc, nc, rhod) >= 0 * si::metres  && "mean droplet radius cannot be < 0");
              assert(tmp >= 0 / si::kilograms / si::seconds && "tmp");

              // so that collisions don't take more n_c than there is
              tmp = std::min(tmp, (nc / si::kilograms) / (dt * si::seconds));
              local_dot_nc -= tmp * si::kilograms * si::seconds;
            }
          }

          dot_rc += local_dot_rc;
          dot_rr += local_dot_rr;
          dot_nc += local_dot_nc;
          dot_nr += local_dot_nr;
        }
      }
    }
  };
};
