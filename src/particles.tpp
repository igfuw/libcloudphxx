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
#include "detail/kernel_utils.hpp"

//kernel definitions
#include "detail/kernel_definitions/hall_efficiencies.hpp"
#include "detail/kernel_definitions/hall_davis_no_waals_efficiencies.hpp"

#include "particles_impl_kernel.ipp"

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
#include "particles_impl_init_grid.ipp"
#include "particles_impl_init_sync.ipp"
#include "particles_impl_init_hskpng.ipp"
#include "particles_impl_init_chem.ipp"
#include "particles_impl_init_kernel.ipp"

#include "particles_impl_hskpng_ijk.ipp"
#include "particles_impl_hskpng_Tpr.ipp"
#include "particles_impl_hskpng_vterm.ipp"
#include "particles_impl_hskpng_sort.ipp"
#include "particles_impl_hskpng_count.ipp"

#include "particles_impl_moms.ipp"
#include "particles_impl_mass_dens.ipp"

#include "particles_impl_fill_outbuf.ipp"

#include "particles_impl_sync.ipp"

#include "particles_impl_adve.ipp"
#include "particles_impl_cond.ipp"
#include "particles_impl_sedi.ipp"
#include "particles_impl_coal.ipp"
#include "particles_impl_chem.ipp"
#include "particles_impl_rcyc.ipp"
#include "particles_impl_bcnd.ipp"
#include "particles_impl_sstp.ipp"
#include "particles_impl_kernel_interpolation.ipp"

