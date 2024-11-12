/** @file
  * @copyright University of Warsaw
  * @brief Rain sedimentation representation for single-moment bulk microphysics
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
    // expects the arguments to be columns with begin() pointing to the lowest level
    // returns rain flux out of the domain
//<listing>
    template <typename real_t, class cont_t>
    real_t rhs_columnwise(
      const opts_t<real_t> &opts,
      cont_t &dot_r_cont,
      const cont_t &rhod_cont,
      const cont_t &r_cont,
      const real_t &dz,
      const std::string& precip_type = "rain"
    )
//</listing>
    {
      using flux_t = quantity<divide_typeof_helper<si::mass_density, si::time>::type, real_t>;

      auto dot_r_unit = si::hertz;

      if (!opts.sedi) return 0;

      //
      flux_t flux_in = 0 * si::kilograms / si::cubic_metres / si::seconds;
      real_t *dot_r = NULL;
      const real_t zero = 0;

      // this should give zero flux from above the domain top
      const real_t *rhod = &*(--(rhod_cont.end())), *r = &zero;

      auto iter = zip(dot_r_cont, rhod_cont, r_cont);
      for (auto tup_ptr = iter.end(); tup_ptr != iter.begin();)
      {
        --tup_ptr;

        const real_t
          *rhod_below  = &std::get<1>(*tup_ptr),
          *r_below     = &std::get<2>(*tup_ptr);

        if (dot_r != NULL) // i.e. all but first (top) grid cell
        {
          flux_t flux_out = real_t(0) * si::kilograms / si::cubic_metres / si::seconds;
          if (precip_type == "rain")
          {
            // terminal momenta at grid-cell edge (to assure precip mass conservation)
            flux_out += -real_t(.5) * ( // averaging + axis orientation
              (*rhod_below * si::kilograms / si::cubic_metres) * formulae::v_term(
                *r_below          * si::kilograms / si::kilograms,
                *rhod_below        * si::kilograms / si::cubic_metres,
                *rhod_cont.begin() * si::kilograms / si::cubic_metres
              ) +
              (*rhod * si::kilograms / si::cubic_metres) * formulae::v_term(
                *r                * si::kilograms / si::kilograms,
                *rhod              * si::kilograms / si::cubic_metres,
                *rhod_cont.begin() * si::kilograms / si::cubic_metres
              )
            ) * (*r * si::kilograms / si::kilograms) / (dz * si::metres);
          }
          else if (precip_type == "iceA")
          {
            // terminal momenta at grid-cell edge (to assure precip mass conservation)
            flux_out += -real_t(.5) * ( // averaging + axis orientation
              (*rhod_below * si::kilograms / si::cubic_metres) * formulae::velocity_iceA(
                *r_below          * si::kilograms / si::kilograms,
                *rhod_below        * si::kilograms / si::cubic_metres
              ) +
              (*rhod * si::kilograms / si::cubic_metres) * formulae::velocity_iceA(
                *r                * si::kilograms / si::kilograms,
                *rhod              * si::kilograms / si::cubic_metres
              )
            ) * (*r * si::kilograms / si::kilograms) / (dz * si::metres);
          }
          else if (precip_type == "iceB")
          {
            // terminal momenta at grid-cell edge (to assure precip mass conservation)
            flux_out += -real_t(.5) * ( // averaging + axis orientation
              (*rhod_below * si::kilograms / si::cubic_metres) * formulae::velocity_iceB(
                *r_below          * si::kilograms / si::kilograms,
                *rhod_below        * si::kilograms / si::cubic_metres
              ) +
              (*rhod * si::kilograms / si::cubic_metres) * formulae::velocity_iceB(
                *r                * si::kilograms / si::kilograms,
                *rhod              * si::kilograms / si::cubic_metres
              )
            ) * (*r * si::kilograms / si::kilograms) / (dz * si::metres);
          }

          *dot_r -= (flux_in - flux_out) / (*rhod * si::kilograms / si::cubic_metres) / dot_r_unit;
          flux_in = flux_out; // inflow = outflow from above
        }

        dot_r = &std::get<0>(*tup_ptr);
        rhod   = rhod_below;
        r     = r_below;
      }

      // the bottom grid cell (with mid-cell vterm approximation)
      flux_t flux_out = real_t(0) * si::kilograms / si::cubic_metres / si::seconds;
      if (precip_type == "rain")
      {
        flux_out += - (*rhod * si::kilograms / si::cubic_metres) * formulae::v_term(
          *r                * si::kilograms / si::kilograms,
          *rhod              * si::kilograms / si::cubic_metres,
          *rhod_cont.begin() * si::kilograms / si::cubic_metres
        ) * (*r * si::kilograms / si::kilograms) / (dz * si::metres);
      }
      else if (precip_type == "iceA")
      {
        flux_out += - (*rhod * si::kilograms / si::cubic_metres) * formulae::velocity_iceA(
          *r                * si::kilograms / si::kilograms,
          *rhod              * si::kilograms / si::cubic_metres
        ) * (*r * si::kilograms / si::kilograms) / (dz * si::metres);
      }
      else if (precip_type == "iceB")
      {
        flux_out += - (*rhod * si::kilograms / si::cubic_metres) * formulae::velocity_iceB(
          *r                * si::kilograms / si::kilograms,
          *rhod              * si::kilograms / si::cubic_metres
        ) * (*r * si::kilograms / si::kilograms) / (dz * si::metres);
      }

      *dot_r -= (flux_in - flux_out) / (*rhod * si::kilograms / si::cubic_metres) / dot_r_unit;
      // outflow from the domain
      return real_t(flux_out / (si::kilograms / si::cubic_metres / si::seconds));
    }
  };
};
