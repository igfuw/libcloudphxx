/** @file
  * @copyright University of Warsaw
  * @brief Rain sedimentation representation for single-moment bulk microphysics
  *   using forcing terms based on the upstrem advection scheme 
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <algorithm>

#include <libcloudph++/blk_1m/formulae.hpp>

namespace libcloudphxx
{
  namespace blk_1m
  {
    // expects the arguments to be columns with begin() pointing to the lowest level
    // returns rain flux out of the domain
//<listing>
    template <typename real_t, class cont_t>
    real_t forcings_columnwise(
      const opts_t<real_t> &opt,
      cont_t &dot_rhod_rr_cont,
      const cont_t &rhod_cont,   
      const cont_t &rhod_rr_cont,
      const real_t &dz // TODO: move to opts
    )   
//</listing>
    {
      using flux_t = quantity<divide_typeof_helper<si::mass_density, si::time>::type, real_t>;

      if (!opt.sedi) return 0;

      // 
      flux_t flux_in = 0 * si::kilograms / si::cubic_metres / si::seconds;
      real_t *dot_rhod_rr = NULL;
      const real_t zero = 0;
      const real_t *rhod, *rhod_rr = &zero;

      auto iter = zip(dot_rhod_rr_cont, rhod_cont, rhod_rr_cont);
      for (auto tup_ptr = iter.end(); tup_ptr != iter.begin();)
      {
        --tup_ptr;

        const real_t
          *rhod_below     = &boost::get<1>(*tup_ptr),
          *rhod_rr_below  = &boost::get<2>(*tup_ptr);

        if (dot_rhod_rr != NULL) // i.e. all but first (top) grid cell
        {
          // terminal velocities at grid-cell edge (to assure precip mass conservation)
	  flux_t flux_out = -.5 * ( // averaging + axis orientation
	    formulae::v_term(
              *rhod_rr_below     * si::kilograms / si::cubic_metres, 
              *rhod_below        * si::kilograms / si::cubic_metres, 
              *rhod_cont.begin() * si::kilograms / si::cubic_metres
            ) + 
	    formulae::v_term(
              *rhod_rr           * si::kilograms / si::cubic_metres,    
              *rhod              * si::kilograms / si::cubic_metres, 
              *rhod_cont.begin() * si::kilograms / si::cubic_metres
            )
	  ) * (*rhod_rr * si::kilograms / si::cubic_metres) / (dz * si::metres);

	  *dot_rhod_rr -= (flux_in - flux_out) * si::seconds * si::cubic_metres / si::kilograms;
          flux_in = flux_out; // inflow = outflow from above
        }

        dot_rhod_rr = &boost::get<0>(*tup_ptr);
         rhod    = rhod_below;
         rhod_rr = rhod_rr_below;
      }
      // the bottom grid cell (with mid-cell vterm approximation)
      flux_t flux_out = - formulae::v_term(
	*rhod_rr           * si::kilograms / si::cubic_metres,    
	*rhod              * si::kilograms / si::cubic_metres, 
	*rhod_cont.begin() * si::kilograms / si::cubic_metres
      ) * (*rhod_rr * si::kilograms / si::cubic_metres) / (dz * si::metres);
      *dot_rhod_rr -= (flux_in - flux_out) * si::seconds * si::cubic_metres / si::kilograms;
      // outflow from the domain
      return real_t(flux_out / (si::kilograms / si::cubic_metres / si::seconds));
    }    
  };
};
