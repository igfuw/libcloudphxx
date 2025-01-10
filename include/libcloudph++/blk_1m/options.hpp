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
      int nwtrph_iters = 3; // number of iterations in Newton-Raphson saturation adjustment
    };
//</listing>
  }; 
};
