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
    template <typename real_t, class container_t>
    quantity<divite_typeof_helper<si::mass_density, si::time>::type, real_t> forcings_columnwise(
      const opts_t<real_t> &opt,
      container_t drhod_rr_cont,
      const container_t rhod_cont,   
      const container_t rhod_rr_cont,
      real_t dz
    )   
    {
      if (!opt.sedi) return;

      // 
      quantity<divide_typeof_helper<si::mass_density, si::time>::type, real_t> 
        flux_in = 0 * si::kilograms / si::cubic_metres / si::seconds;
      real_t *drhod_rr = NULL;
      const real_t zero = 0;
      const real_t *rhod, *rhod_rr = &zero;

      auto iter = zip(drhod_rr_cont, rhod_cont, rhod_rr_cont);
      for (auto tup_ptr = --iter.end(); tup_ptr != --iter.begin(); --tup_ptr)
      {
        const real_t
          *rhod_below     = &boost::get<1>(*tup_ptr),
          *rhod_rr_below  = &boost::get<2>(*tup_ptr);

        if (drhod_rr != NULL) // i.e. all but first (top) grid cell
        {
          // terminal velocities at grid-cell edge (to assure precip mass conservation)
	  quantity<divide_typeof_helper<si::mass_density, si::time>::type, real_t> flux_out = -.5 * ( // averaging + axis orientation
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

	  *drhod_rr -= (flux_in - flux_out) * si::seconds * si::cubic_metres / si::kilograms;
          flux_in = flux_out; // inflow = outflow from above
        }

        drhod_rr = &boost::get<0>(*tup_ptr);
         rhod    = rhod_below;
         rhod_rr = rhod_rr_below;
      }
      // assumption: inflow to the bottom grid cell = outflow from the domain
      return flux_in; // (was: *drhod_rr -= ...)
    }    
  }
};
