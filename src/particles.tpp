// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief Thrust-based CPU/GPU particle-tracking logic for Lagrangian microphysics
  */

#include <iostream>

#include <libcloudph++/lgrngn/particles.hpp>

#include "detail/config.hpp"
#include "detail/thrust.hpp"
#include "detail/urand.hpp"
#include "detail/eval_and_oper.hpp"
#include "detail/kernel_utils.hpp"
#include "detail/wang_collision_enhancement.hpp"
#include "detail/kernel_onishi_nograv.hpp"
#include "detail/checknan.hpp"
#include "detail/formatter.cpp"

//kernel definitions
#include "detail/kernel_definitions/hall_efficiencies.hpp"
#include "detail/kernel_definitions/hall_davis_no_waals_efficiencies.hpp"
#include "detail/kernel_definitions/vohl_davis_no_waals_efficiencies.hpp"
#include "detail/kernel_definitions/hall_pinsky_stratocumulus_efficiencies.hpp"
#include "detail/kernel_definitions/hall_pinsky_cumulonimbus_efficiencies.hpp"
#include "detail/kernel_definitions/hall_pinsky_1000mb_grav_efficiencies.hpp"

#include "particles_impl_kernel.ipp"

// public API
#include "particles_pimpl_ctor.ipp"
#include "particles_init.ipp"
#include "particles_step.ipp"
#include "particles_diag.ipp"

// details
#include "particles_impl_dist_analysis.ipp"
#include "particles_impl_init_SD_with_distros_sd_conc.ipp"
#include "particles_impl_init_SD_with_distros_tail.ipp"
#include "particles_impl_init_SD_with_distros_const_multi.ipp"
#include "particles_impl_init_SD_with_distros.ipp"
#include "particles_impl_init_SD_with_sizes.ipp"
#include "particles_impl_init_dry_sd_conc.ipp"
#include "particles_impl_init_dry_const_multi.ipp"
#include "particles_impl_init_dry_dry_sizes.ipp"
#include "particles_impl_init_kappa.ipp"
#include "particles_impl_init_n.ipp"
#include "particles_impl_init_wet.ipp"
#include "particles_impl_init_xyz.ipp"
#include "particles_impl_init_ijk.ipp"
#include "particles_impl_init_count_num.ipp"
#include "particles_impl_init_e2l.ipp"
#include "particles_impl_init_grid.ipp"
#include "particles_impl_init_sync.ipp"
#include "particles_impl_init_hskpng_npart.ipp"
#include "particles_impl_init_hskpng_ncell.ipp"
#include "particles_impl_init_chem.ipp"
#include "particles_impl_init_kernel.ipp"
#include "particles_impl_step_finalize.ipp"
#include "particles_impl_init_vterm.ipp"
#include "particles_impl_init_sanity_check.ipp"
#include "particles_impl_update_th_rv.ipp"

#include "particles_impl_hskpng_ijk.ipp"
#include "particles_impl_hskpng_Tpr.ipp"
#include "particles_impl_hskpng_vterm.ipp"
#include "particles_impl_hskpng_sort.ipp"
#include "particles_impl_hskpng_count.ipp"
#include "particles_impl_hskpng_remove.ipp"
#include "particles_impl_hskpng_resize.ipp"

#include "particles_impl_moms.ipp"
#include "particles_impl_mass_dens.ipp"

#include "particles_impl_fill_outbuf.ipp"

#include "particles_impl_sync.ipp"

#include "particles_impl_bcnd.ipp" // bcnd has to be b4 adve for periodic struct; move it to separate file in detail...
#include "particles_impl_adve.ipp"
#include "particles_impl_cond_common.ipp"
#include "particles_impl_cond.ipp"
#include "particles_impl_cond_sstp.ipp"
#include "particles_impl_sedi.ipp"
#include "particles_impl_coal.ipp"
#include "particles_impl_chem_ante.ipp"
#include "particles_impl_chem_henry.ipp"
#include "particles_impl_chem_dissoc.ipp"
#include "particles_impl_chem_strength.ipp"
#include "particles_impl_chem_react.ipp"
#include "particles_impl_rcyc.ipp"
#include "particles_impl_sstp.ipp"
#include "particles_impl_sstp_chem.ipp"
#include "particles_impl_src.ipp"
#include "particles_impl_kernel_interpolation.ipp"

