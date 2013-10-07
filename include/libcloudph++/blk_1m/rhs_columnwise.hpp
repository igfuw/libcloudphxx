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
    real_t rhs_columnwise(
      const opts_t<real_t> &opts,
      cont_t &dot_rho_r_cont,
      const cont_t &rho_d_cont,   
      const cont_t &rho_r_cont,
      const real_t &dz 
    )   
//</listing>
    {
      using flux_t = quantity<divide_typeof_helper<si::mass_density, si::time>::type, real_t>;

      if (!opts.sedi) return 0;

      // 
      flux_t flux_in = 0 * si::kilograms / si::cubic_metres / si::seconds;
      real_t *dot_rho_r = NULL;
      const real_t zero = 0;
  
      // this should give zero flux from above the domain top
      const real_t *rho_d = &*(--(rho_d_cont.end())), *rho_r = &zero;

      auto iter = zip(dot_rho_r_cont, rho_d_cont, rho_r_cont);
      for (auto tup_ptr = iter.end(); tup_ptr != iter.begin();)
      {
        --tup_ptr;

        const real_t
          *rho_d_below  = &boost::get<1>(*tup_ptr),
          *rho_r_below  = &boost::get<2>(*tup_ptr);

        if (dot_rho_r != NULL) // i.e. all but first (top) grid cell
        {
          // terminal velocities at grid-cell edge (to assure precip mass conservation)
	  flux_t flux_out = -real_t(.5) * ( // averaging + axis orientation
	    formulae::v_term(
              *rho_r_below        * si::kilograms / si::cubic_metres, 
              *rho_d_below        * si::kilograms / si::cubic_metres, 
              *rho_d_cont.begin() * si::kilograms / si::cubic_metres
            ) + 
	    formulae::v_term(
              *rho_r              * si::kilograms / si::cubic_metres,    
              *rho_d              * si::kilograms / si::cubic_metres, 
              *rho_d_cont.begin() * si::kilograms / si::cubic_metres
            )
	  ) * (*rho_r * si::kilograms / si::cubic_metres) / (dz * si::metres);

	  *dot_rho_r -= (flux_in - flux_out) * si::seconds * si::cubic_metres / si::kilograms;
          flux_in = flux_out; // inflow = outflow from above
        }

        dot_rho_r = &boost::get<0>(*tup_ptr);
        rho_d     = rho_d_below;
        rho_r     = rho_r_below;
      }
      // the bottom grid cell (with mid-cell vterm approximation)
      flux_t flux_out = - formulae::v_term(
	*rho_r              * si::kilograms / si::cubic_metres,    
	*rho_d              * si::kilograms / si::cubic_metres, 
	*rho_d_cont.begin() * si::kilograms / si::cubic_metres
      ) * (*rho_r * si::kilograms / si::cubic_metres) / (dz * si::metres);
      *dot_rho_r -= (flux_in - flux_out) * si::seconds * si::cubic_metres / si::kilograms;
      // outflow from the domain
      return real_t(flux_out / (si::kilograms / si::cubic_metres / si::seconds));
    }    
  };
};
