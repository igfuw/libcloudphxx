/** @file
  * @copyright University of Warsaw
  * @brief Rain and ice sedimentation representation for single-moment bulk microphysics
  *   using forcing terms based on the upstrem advection scheme
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include "extincl.hpp"

namespace libcloudphxx
{
  namespace blk_1m
  {
    enum class ice_t {iceA, iceB};

    // expects the arguments to be columns with begin() pointing to the lowest level
    // returns rain flux out of the domain
//<listing>
    template <typename real_t, class cont_t>
    real_t rhs_columnwise(
      const opts_t<real_t> &opts,
      cont_t &dot_rr_cont,
      const cont_t &rhod_cont,
      const cont_t &rr_cont,
      const real_t &dz
    )
//</listing>
    {
      using flux_t = quantity<divide_typeof_helper<si::mass_density, si::time>::type, real_t>;

      auto dot_rr_unit = si::hertz;

      if (!opts.sedi) return 0;

      //
      flux_t flux_in = 0 * si::kilograms / si::cubic_metres / si::seconds;
      real_t *dot_rr = NULL;
      const real_t zero = 0;

      // this should give zero flux from above the domain top
      const real_t *rhod = &*(--(rhod_cont.end())), *rr = &zero;

      auto iter = zip(dot_rr_cont, rhod_cont, rr_cont);
      for (auto tup_ptr = iter.end(); tup_ptr != iter.begin();)
      {
        --tup_ptr;

        const real_t
          *rhod_below  = &std::get<1>(*tup_ptr),
          *rr_below    = &std::get<2>(*tup_ptr);

        if (dot_rr != NULL) // i.e. all but first (top) grid cell
        {
          // terminal momenta at grid-cell edge (to assure precip mass conservation)
          flux_t flux_out = -real_t(.5) * ( // averaging + axis orientation
            (*rhod_below * si::kilograms / si::cubic_metres) * formulae::v_term(
              *rr_below          * si::kilograms / si::kilograms,
              *rhod_below        * si::kilograms / si::cubic_metres,
              *rhod_cont.begin() * si::kilograms / si::cubic_metres
            ) +
            (*rhod * si::kilograms / si::cubic_metres) * formulae::v_term(
              *rr                * si::kilograms / si::kilograms,
              *rhod              * si::kilograms / si::cubic_metres,
              *rhod_cont.begin() * si::kilograms / si::cubic_metres
            )
          ) * (*rr * si::kilograms / si::kilograms) / (dz * si::metres);

          *dot_rr -= (flux_in - flux_out) / (*rhod * si::kilograms / si::cubic_metres) / dot_rr_unit;
          flux_in = flux_out; // inflow = outflow from above
        }

        dot_rr = &std::get<0>(*tup_ptr);
        rhod   = rhod_below;
        rr     = rr_below;
      }

      // the bottom grid cell (with mid-cell vterm approximation)
      flux_t flux_out = - (*rhod * si::kilograms / si::cubic_metres) * formulae::v_term(
        *rr                * si::kilograms / si::kilograms,
        *rhod              * si::kilograms / si::cubic_metres,
        *rhod_cont.begin() * si::kilograms / si::cubic_metres
      ) * (*rr * si::kilograms / si::kilograms) / (dz * si::metres);
      *dot_rr -= (flux_in - flux_out) / (*rhod * si::kilograms / si::cubic_metres) / dot_rr_unit;

      // outflow from the domain
      return real_t(flux_out / (si::kilograms / si::cubic_metres / si::seconds));
    }


//<listing>
    template <typename real_t, class cont_t>
    real_t rhs_columnwise_ice(
      const opts_t<real_t> &opts,
      cont_t &dot_ri_cont,
      const cont_t &rhod_cont,
      const cont_t &ri_cont,
      const real_t &dz,
      const ice_t& ice_type
    )
//</listing>
    {
      using flux_t = quantity<divide_typeof_helper<si::mass_density, si::time>::type, real_t>;

      auto dot_ri_unit = si::hertz;

      if (!opts.sedi) return 0;

      //
      flux_t flux_in = 0 * si::kilograms / si::cubic_metres / si::seconds;
      real_t *dot_ri = NULL;
      const real_t zero = 0;

      // this should give zero flux from above the domain top
      const real_t *rhod = &*(--(rhod_cont.end())), *ri = &zero;

      auto iter = zip(dot_ri_cont, rhod_cont, ri_cont);
      for (auto tup_ptr = iter.end(); tup_ptr != iter.begin();)
      {
        --tup_ptr;

        const real_t
          *rhod_below  = &std::get<1>(*tup_ptr),
          *ri_below     = &std::get<2>(*tup_ptr);

        if (dot_ri != NULL) // i.e. all but first (top) grid cell
        {
          flux_t flux_out = real_t(0) * si::kilograms / si::cubic_metres / si::seconds;
          if (ice_type == ice_t::iceA)
          {
            // terminal momenta at grid-cell edge (to assure precip mass conservation)
            flux_out += -real_t(.5) * ( // averaging + axis orientation
              (*rhod_below * si::kilograms / si::cubic_metres) * formulae::velocity_iceA(
                *ri_below          * si::kilograms / si::kilograms,
                *rhod_below        * si::kilograms / si::cubic_metres
              ) +
              (*rhod * si::kilograms / si::cubic_metres) * formulae::velocity_iceA(
                *ri                * si::kilograms / si::kilograms,
                *rhod              * si::kilograms / si::cubic_metres
              )
            ) * (*ri * si::kilograms / si::kilograms) / (dz * si::metres);
          }
          else if (ice_type == ice_t::iceB)
          {
            // terminal momenta at grid-cell edge (to assure precip mass conservation)
            flux_out += -real_t(.5) * ( // averaging + axis orientation
              (*rhod_below * si::kilograms / si::cubic_metres) * formulae::velocity_iceB(
                *ri_below          * si::kilograms / si::kilograms,
                *rhod_below        * si::kilograms / si::cubic_metres
              ) +
              (*rhod * si::kilograms / si::cubic_metres) * formulae::velocity_iceB(
                *ri                * si::kilograms / si::kilograms,
                *rhod              * si::kilograms / si::cubic_metres
              )
            ) * (*ri * si::kilograms / si::kilograms) / (dz * si::metres);
          }

          *dot_ri -= (flux_in - flux_out) / (*rhod * si::kilograms / si::cubic_metres) / dot_ri_unit;
          flux_in = flux_out; // inflow = outflow from above
        }

        dot_ri = &std::get<0>(*tup_ptr);
        rhod   = rhod_below;
        ri     = ri_below;
      }

      // the bottom grid cell (with mid-cell vterm approximation)
      flux_t flux_out = real_t(0) * si::kilograms / si::cubic_metres / si::seconds;
      if (ice_type == ice_t::iceA)
      {
        flux_out += - (*rhod * si::kilograms / si::cubic_metres) * formulae::velocity_iceA(
          *ri                * si::kilograms / si::kilograms,
          *rhod              * si::kilograms / si::cubic_metres
        ) * (*ri * si::kilograms / si::kilograms) / (dz * si::metres);
      }
      else if (ice_type == ice_t::iceB)
      {
        flux_out += - (*rhod * si::kilograms / si::cubic_metres) * formulae::velocity_iceB(
          *ri                * si::kilograms / si::kilograms,
          *rhod              * si::kilograms / si::cubic_metres
        ) * (*ri * si::kilograms / si::kilograms) / (dz * si::metres);
      }

      *dot_ri -= (flux_in - flux_out) / (*rhod * si::kilograms / si::cubic_metres) / dot_ri_unit;
      // outflow from the domain
      return real_t(flux_out / (si::kilograms / si::cubic_metres / si::seconds));
    }
  };
};
