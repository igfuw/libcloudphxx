// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief Thrust-based CPU/GPU particle-tracking logic for Lagrangian microphysics
  */

#include <iostream>

#include "../include/lgrngn/particles.hpp"

#include "detail/urand.hpp"
#include "detail/thrust.hpp"

#include "particles_pimpl_ctor.ipp"
#include "particles_init.ipp"
