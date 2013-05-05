#pragma once
#include <libcloudph++/common/phc_kessler.hpp>
#include <algorithm>

namespace libcloudphxx
{
  namespace bulk
  {
    // expects the arguments to be columns with begin() pointing to the lowest level
    template <typename real_t, class container_t>
    void sediment(
      container_t drhod_rr_cont,
      const container_t rhod_cont,   
      const container_t rhod_rr_cont,
      real_t dz
    )   
    {
      // TODO: return accumulated rainfall?

      // 
      quantity<divide_typeof_helper<si::mass_density, si::time>::type, real_t> flux_in = 0 * si::kilograms / si::cubic_metres / si::seconds;
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
	    phc::v_term(
              *rhod_rr_below     * si::kilograms / si::cubic_metres, 
              *rhod_below        * si::kilograms / si::cubic_metres, 
              *rhod_cont.begin() * si::kilograms / si::cubic_metres
            ) + 
	    phc::v_term(
              *rhod_rr           * si::kilograms / si::cubic_metres,    
              *rhod              * si::kilograms / si::cubic_metres, 
              *rhod_cont.begin() * si::kilograms / si::cubic_metres
            )
	  ) * (*rhod_rr * si::kilograms / si::cubic_metres) / (dz * si::metres);

	  *drhod_rr -= (flux_in - flux_out) * si::seconds * si::cubic_metres / si::kilograms;
          flux_in = flux_out; // inflow = outflow from above
        }

        drhod_rr = &boost::get<0>(*tup_ptr);
        rhod = rhod_below;
        rhod_rr = rhod_rr_below;
      }
      // inflow to the bottom grid cell
      *drhod_rr -= flux_in * si::seconds * si::cubic_metres / si::kilograms;
      // TODO: outflow from the bottom grid-cell
    }    
  }
};
