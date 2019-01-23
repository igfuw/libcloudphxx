/** @file
  * @copyright University of Warsaw
  * @brief Rain sedimentation representation for single-moment bulk microphysics
  *   using forcing terms based on the upstrem advection scheme
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <libcloudph++/blk_2m/extincl.hpp>

namespace libcloudphxx
{
  namespace blk_2m
  {
    // expects the arguments to be columns with begin() pointing to the lowest level
    // returns rain flux out of the domain
//<listing>
    template <typename real_t, class cont_t>
    real_t rhs_columnwise(
      const opts_t<real_t> &opts,
      cont_t &dot_rr_cont,
      cont_t &dot_nr_cont,
      const cont_t &rhod_cont,
      const cont_t &rr_cont,
      const cont_t &nr_cont,
      const real_t &dt,
      const real_t &dz
    )
//</listing>
    {
      if (!opts.sedi) return 0;

      using flux_rr = quantity<divide_typeof_helper<si::mass_density, si::time>::type, real_t>;
      using flux_nr = quantity<divide_typeof_helper<si::frequency, si::volume>::type, real_t>;

      auto dot_rr_unit = si::hertz;
      auto dot_nr_unit = si::hertz / si::kilograms;

      auto rflux_unit = si::kilograms / si::seconds / si::cubic_metres;
      auto nflux_unit = si::hertz / si::cubic_metres;

      flux_rr flux_rr_in = 0 * si::kilograms / si::cubic_metres / si::seconds;
      flux_nr flux_nr_in = 0 / si::cubic_metres / si::seconds;

      real_t *dot_rr = NULL;
      real_t *dot_nr = NULL;
      const real_t zero = 0;

      // initial values that should give zero flux from above the domain top
      const real_t
        *rr = &zero,
        *nr = &zero,
        *rhod = &*(--(rhod_cont.end()));

      auto iter = zip(rhod_cont, rr_cont, nr_cont, dot_rr_cont, dot_nr_cont);
      for (auto tup_ptr = iter.end(); tup_ptr != iter.begin();)
      {
        --tup_ptr;

        const real_t
          *rhod_below = &std::get<0>(*tup_ptr),
          *rr_below   = &std::get<1>(*tup_ptr),
          *nr_below   = &std::get<2>(*tup_ptr);

        if (dot_rr != NULL) // i.e. all but first (top) grid cell
        {
          // terminal velocities at grid-cell edge (to assure precip mass conservation)
          quantity<multiply_typeof_helper<si::velocity, si::mass_density>::type, real_t> tmp_mom_m  = -real_t(.5) * ( // averaging + axis orientation
	    (*rhod_below * si::kilograms / si::cubic_metres) * formulae::v_term_m(
              *rhod_below * si::kilograms / si::cubic_metres,
              *rr_below * si::kilograms / si::kilograms,
              *nr_below / si::kilograms
            ) +
	    (*rhod * si::kilograms / si::cubic_metres) * formulae::v_term_m(
              *rhod * si::kilograms / si::cubic_metres,
              *rr * si::kilograms / si::kilograms,
              *nr / si::kilograms
            )
	  );

          quantity<multiply_typeof_helper<si::velocity, si::mass_density>::type, real_t> tmp_mom_n  = -real_t(.5) * ( // averaging + axis orientation
	    (*rhod_below * si::kilograms / si::cubic_metres) * formulae::v_term_n(
              *rhod_below * si::kilograms / si::cubic_metres,
              *rr_below * si::kilograms / si::kilograms,
              *nr_below / si::kilograms
            ) +
	    (*rhod * si::kilograms / si::cubic_metres) * formulae::v_term_n(
              *rhod * si::kilograms / si::cubic_metres,
              *rr * si::kilograms / si::kilograms,
              *nr / si::kilograms
            )
	  );

          flux_rr flux_rr_out = tmp_mom_m * (*rr * si::kilograms / si::kilograms) / (dz * si::metres);
          flux_rr_out = - std::min(real_t(-flux_rr_out / rflux_unit), *rhod * (*rr + dt * *dot_rr) / dt) * rflux_unit;

          flux_nr flux_nr_out = tmp_mom_n * (*nr / si::kilograms) / (dz * si::metres);
          flux_nr_out = - std::min(real_t(-flux_nr_out / nflux_unit), *rhod * (*nr + dt * *dot_nr) / dt) * nflux_unit;

	  *dot_rr -= (flux_rr_in - flux_rr_out) / (*rhod * si::kilograms / si::cubic_metres) / dot_rr_unit;
          flux_rr_in = flux_rr_out; // inflow = outflow from above
	  *dot_nr -= (flux_nr_in - flux_nr_out) / (*rhod * si::kilograms / si::cubic_metres) / dot_nr_unit;
          flux_nr_in = flux_nr_out; // inflow = outflow from above
        }

        dot_rr = &std::get<3>(*tup_ptr);
        dot_nr = &std::get<4>(*tup_ptr);
        rhod = rhod_below;
        rr = rr_below;
        nr = nr_below;
      }

      // the bottom grid cell (with mid-cell vterm approximation)
      quantity<multiply_typeof_helper<si::velocity, si::mass_density>::type, real_t>
        tmp_mom_m = - (*rhod * si::kilograms / si::cubic_metres) * formulae::v_term_m(
	  *rhod * si::kilograms / si::cubic_metres,
	  *rr * si::kilograms / si::kilograms,
	  *nr / si::kilograms
	);
      quantity<multiply_typeof_helper<si::velocity, si::mass_density>::type, real_t>
	tmp_mom_n = - (*rhod * si::kilograms / si::cubic_metres) * formulae::v_term_n(
	  *rhod * si::kilograms / si::cubic_metres,
	  *rr * si::kilograms / si::kilograms,
	  *nr / si::kilograms
	);

      // outflow from the domain
      {
        flux_nr flux_nr_out = tmp_mom_n * (*nr / si::kilograms) / (dz * si::metres);
        flux_nr_out = - std::min(real_t(-flux_nr_out / nflux_unit), *rhod * (*nr + dt * *dot_nr) / dt) * nflux_unit;
        *dot_nr -= (flux_nr_in - flux_nr_out) / (*rhod * si::kilograms / si::cubic_metres) / dot_nr_unit;
      }

      {
        flux_rr flux_rr_out = tmp_mom_m * (*rr * si::kilograms / si::kilograms) / (dz * si::metres);
        flux_rr_out = - std::min(real_t(-flux_rr_out / rflux_unit), *rhod * (*rr + dt * *dot_rr) / dt) * rflux_unit;
        *dot_rr -= (flux_rr_in - flux_rr_out) / (*rhod * si::kilograms / si::cubic_metres) / dot_rr_unit;

        return flux_rr_out / (si::kilograms / si::cubic_metres / si::seconds);
      }
    }
  };
};
