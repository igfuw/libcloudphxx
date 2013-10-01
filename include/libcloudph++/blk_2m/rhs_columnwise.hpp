/** @file
  * @copyright University of Warsaw
  * @brief Rain sedimentation representation for single-moment bulk microphysics
  *   using forcing terms based on the upstrem advection scheme 
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

#include <algorithm>

#include <libcloudph++/common/detail/zip.hpp>
#include <libcloudph++/blk_2m/terminal_vel_formulae.hpp>  //TODO so far the same parametrisation as in blk_1m 

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
      cont_t &dot_rho_r_cont,
      cont_t &dot_n_r_cont,
      const cont_t &rho_r_cont,
      const cont_t &n_r_cont,
      const real_t &dt,
      const real_t &dz
    )   
//</listing>
    {
      using flux_rr = quantity<divide_typeof_helper<si::mass_density, si::time>::type, real_t>;

      if (!opts.sedi) return 0;

      using flux_nr = quantity<divide_typeof_helper<si::frequency, si::volume>::type, real_t>;
     
      flux_rr flux_rr_in = 0 * si::kilograms / si::cubic_metres / si::seconds;
      flux_nr flux_nr_in = 0 / si::cubic_metres / si::seconds;


      real_t *dot_rho_r = NULL;
      real_t *dot_n_r = NULL;
      const real_t zero = 0;
      const real_t *rho_r, *n_r = &zero;

      auto iter = zip(rho_r_cont, n_r_cont, dot_rho_r_cont, dot_n_r_cont);
      for (auto tup_ptr = iter.end(); tup_ptr != iter.begin();)
      {
        --tup_ptr;

        const real_t
          *rho_r_below  = &boost::get<0>(*tup_ptr),
          *n_r_below  = &boost::get<1>(*tup_ptr);

        if (dot_rho_r != NULL) // i.e. all but first (top) grid cell
        {
          // terminal velocities at grid-cell edge (to assure precip mass conservation)
          quantity<si::velocity, real_t> tmp_vel  = -real_t(.5) * ( // averaging + axis orientation
	    formulae::v_term(
              *rho_r_below * si::kilograms / si::cubic_metres, 
              *n_r_below / si::cubic_metres 
            ) + 
	    formulae::v_term(
              *rho_r * si::kilograms / si::cubic_metres,
              *n_r / si::cubic_metres
            )
	  ); 
          
          auto rflux_unit = si::kilograms / si::seconds / si::cubic_metres;
          auto nflux_unit = si::hertz / si::cubic_metres;

          flux_rr flux_rr_out = tmp_vel * (*rho_r * si::kilograms / si::cubic_metres) / (dz * si::metres);
          flux_rr_out = std::min(real_t(flux_rr_out / rflux_unit), (*rho_r + dt * *dot_rho_r) / dt) * rflux_unit;

          flux_nr flux_nr_out = tmp_vel * (*n_r / si::cubic_metres) / (dz * si::metres);
          flux_nr_out = std::min(real_t(flux_nr_out / nflux_unit), (*n_r + dt * *dot_n_r) / dt) * nflux_unit;

	  *dot_rho_r -= (flux_rr_in - flux_rr_out) / rflux_unit;
          flux_rr_in = flux_rr_out; // inflow = outflow from above
	  *dot_n_r -= (flux_nr_in - flux_nr_out) / nflux_unit;
          flux_nr_in = flux_nr_out; // inflow = outflow from above
        }

        dot_rho_r = &boost::get<2>(*tup_ptr);
        dot_n_r = &boost::get<3>(*tup_ptr);
        rho_r = rho_r_below;
        n_r = n_r_below;
      }

      // the bottom grid cell (with mid-cell vterm approximation)
      quantity<si::velocity, real_t> tmp_vel = - formulae::v_term(
	*rho_r * si::kilograms / si::cubic_metres,
	*n_r / si::cubic_metres
      ); 

      {
        flux_nr flux_nr_out = tmp_vel * (*n_r / si::cubic_metres) / (dz * si::metres);
        *dot_n_r -= (flux_nr_in - flux_nr_out) * si::seconds * si::cubic_metres;
      }

      // outflow from the domain
      {
        flux_rr flux_rr_out = tmp_vel * (*rho_r * si::kilograms / si::cubic_metres) / (dz * si::metres);
        *dot_rho_r -= (flux_rr_in - flux_rr_out) * si::seconds * si::cubic_metres / si::kilograms;

        return flux_rr_out / (si::kilograms / si::cubic_metres / si::seconds);
      }
    }    
  };
};
