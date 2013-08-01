// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief Thrust-based CPU/GPU particle-tracking logic for Lagrangian microphysics
  */

#include <iostream>

#include <libcloudph++/lgrngn/particles.hpp>

#include "detail/thrust.hpp"
#include "detail/urand.hpp"

// public API
#include "particles_pimpl_ctor.ipp"
#include "particles_init.ipp"
#include "particles_step.ipp"
#include "particles_diag.ipp"

// details
#include "particles_impl_init_dry.ipp"
#include "particles_impl_init_wet.ipp"
#include "particles_impl_init_xyz.ipp"
#include "particles_impl_init_e2l.ipp"
#include "particles_impl_init_sync.ipp"
#include "particles_impl_init_hskpng.ipp"

#include "particles_impl_hskpng_ijk.ipp"
#include "particles_impl_hskpng_Tpr.ipp"
#include "particles_impl_hskpng_sort.ipp"
#include "particles_impl_hskpng_count.ipp"

#include "particles_impl_sync.ipp"

#include "particles_impl_adve.ipp"
