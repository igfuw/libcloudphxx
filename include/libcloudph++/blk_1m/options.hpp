/** @file
  * @copyright University of Warsaw
  * @brief Definition of a structure holding options for single-moment bulk microphysics
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

namespace libcloudphxx
{
  namespace blk_1m 
  { 
//<listing>
    template<typename real_t>
    struct opts_t {
      bool 
        cond = true,    // condensation
        cevp = true,    // evaporation of cloud
        revp = true,    // evaporation of rain
        conv = true,    // autoconversion
        accr = true,    // accretion
        sedi = true,    // sedimentation
        homA1 = true,   // homogeneous nucleation of ice A from water vapor
        homA2 = true,   // homogeneous nucleation of ice A from cloud droplets
        hetA = true,    // heterogeneous nucleation of ice A
        hetB = true,    // heterogeneous nucleation of ice B
        depA = true,    // depositional growth of ice A
        depB = true,    // depositional growth of ice B
        rimA = true,    // growth of ice A by riming
        rimB = true,    // growth of ice B by riming
        melA = true,    // melting of ice A
        melB = true;    // melting of ice B
      real_t 
        r_c0   = 5e-4,   // autoconv. threshold
        k_acnv = 0.001,  // Kessler autoconversion (eq. 5a in Grabowski & Smolarkiewicz 1996)
        r_eps  = 2e-5;   // absolute tolerance

      bool adj_nwtrph = true; // if true, use simpler Newton-Raphson iteration in saturation adjustment; otherwise use RK4 from boost.odeint
      int nwtrph_iters = 3; // number of iterations in Newton-Raphson saturation adjustment

      // NOTE: only tested combinations are: th_dry == true && const_p == false; th_dry == false && const_p == true
      // NOTE:  th_dry == true && const_p == false doesn't work very well with Newton-Raphson, e.g. sat_adj_blk_1m test (TODO: probably Newton-Raphson needs to be fixed)
      bool th_dry  = true, // if true, input and output theta are dry-air potential temperature; if false, they are 'standard' potential temperature
           const_p = false; // if true, pressure is equal to a supplied profile except for solving velocity (e.g. anelastic model); if false, pressure comes from the gas equation
    };
//</listing>
  }; 
};
