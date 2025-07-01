// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief Thrust-based CPU/GPU particle-tracking logic for Lagrangian microphysics
  */

#include <iostream>
#if !defined(NDEBUG)
  #include <iomanip>
#endif

#include <libcloudph++/lgrngn/particles.hpp>

#include "detail/config.hpp"
#include "detail/thrust.hpp"
#include "detail/urand.hpp"
#include "detail/eval_and_oper.hpp"
#include "detail/kernel_utils.hpp"
#include "detail/wang_collision_enhancement.hpp"
#include "detail/kernel_onishi_nograv.hpp"
#include "detail/bcond.hpp"
#include "detail/checknan.hpp"
#include "detail/formatter.cpp"
#include "detail/tpl_calc_wrapper.hpp"
#include "detail/kernels.hpp"
#include "detail/kernel_interpolation.hpp"
#include "detail/functors_host.hpp"
#include "detail/ran_with_mpi.hpp"

//kernel definitions
#include "detail/kernel_definitions/hall_efficiencies.hpp"
#include "detail/kernel_definitions/hall_davis_no_waals_efficiencies.hpp"
#include "detail/kernel_definitions/vohl_davis_no_waals_efficiencies.hpp"
#include "detail/kernel_definitions/hall_pinsky_stratocumulus_efficiencies.hpp"
#include "detail/kernel_definitions/hall_pinsky_cumulonimbus_efficiencies.hpp"
#include "detail/kernel_definitions/hall_pinsky_1000mb_grav_efficiencies.hpp"

#if defined(USE_MPI)
  #include <mpi.h>
  // MPI init
  #include "detail/mpi_init.hpp"
  #include "detail/get_mpi_type.hpp"
#endif

// public API
#include "particles_ctor.ipp"
#include "particles_init.ipp"
#include "particles_step.ipp"
#include "particles_diag.ipp"

// details
#include "impl/particles_impl.ipp"
#include "impl/particles_impl_dist_analysis.ipp"
#include "impl/particles_impl_reserve_hskpng_npart.ipp"
#include "impl/particles_impl_init_SD_with_distros_sd_conc.ipp"
#include "impl/particles_impl_init_SD_with_distros_tail.ipp"
#include "impl/particles_impl_init_SD_with_distros_const_multi.ipp"
#include "impl/particles_impl_init_SD_with_distros.ipp"
#include "impl/particles_impl_init_SD_with_sizes.ipp"
#include "impl/particles_impl_init_dry_sd_conc.ipp"
#include "impl/particles_impl_init_dry_const_multi.ipp"
#include "impl/particles_impl_init_dry_dry_sizes.ipp"
#include "impl/particles_impl_init_kappa.ipp"
#include "impl/particles_impl_init_rd3_insol.ipp"
#include "impl/particles_impl_init_T_freeze.ipp"
#include "impl/particles_impl_init_incloud_time.ipp"
#include "impl/particles_impl_init_n.ipp"
#include "impl/particles_impl_init_wet.ipp"
#include "impl/particles_impl_init_xyz.ipp"
#include "impl/particles_impl_init_ijk.ipp"
#include "impl/particles_impl_init_count_num.ipp"
#include "impl/particles_impl_init_e2l.ipp"
#include "impl/particles_impl_init_grid.ipp"
#include "impl/particles_impl_init_sync.ipp"
#include "impl/particles_impl_init_hskpng_ncell.ipp"
#include "impl/particles_impl_init_chem.ipp"
#include "impl/particles_impl_init_kernel.ipp"
#include "impl/particles_impl_post_copy.ipp"
#include "impl/particles_impl_init_vterm.ipp"
#include "impl/particles_impl_init_sanity_check.ipp"
#include "impl/particles_impl_xchng_domains.ipp"
#include "impl/particles_impl_xchng_courants.ipp"
#include "impl/particles_impl_update_th_rv.ipp"
#include "impl/particles_impl_hskpng_ijk.ipp"
#include "impl/particles_impl_hskpng_Tpr.ipp"
#include "impl/particles_impl_hskpng_mfp.ipp"
#include "impl/particles_impl_hskpng_vterm.ipp"
#include "impl/particles_impl_hskpng_turb_vel.ipp"
#include "impl/particles_impl_hskpng_turb_ss.ipp"
#include "impl/particles_impl_hskpng_tke.ipp"
#include "impl/particles_impl_hskpng_sort.ipp"
#include "impl/particles_impl_hskpng_count.ipp"
#include "impl/particles_impl_hskpng_remove.ipp"
#include "impl/particles_impl_hskpng_resize.ipp"
#include "impl/particles_impl_moms.ipp"
#include "impl/particles_impl_mass_dens.ipp"
#include "impl/particles_impl_fill_outbuf.ipp"
#include "impl/particles_impl_distmem_access.ipp"
#include "impl/particles_impl_sync.ipp"
#include "impl/particles_impl_bcnd.ipp" // bcnd has to be b4 adve for periodic struct; move it to separate file in detail...
#include "impl/particles_impl_adve.ipp"
#include "impl/particles_impl_turb_adve.ipp"
#include "impl/particles_impl_cond_common.ipp"
#include "impl/particles_impl_cond.ipp"
#include "impl/particles_impl_cond_sstp.ipp"
#include "impl/particles_impl_nucleation.ipp"
#include "impl/particles_impl_sedi.ipp"
#include "impl/particles_impl_subs.ipp"
#include "impl/particles_impl_coal.ipp"
#include "impl/particles_impl_chem_ante.ipp"
#include "impl/particles_impl_chem_henry.ipp"
#include "impl/particles_impl_chem_dissoc.ipp"
#include "impl/particles_impl_chem_strength.ipp"
#include "impl/particles_impl_chem_react.ipp"
#include "impl/particles_impl_rcyc.ipp"
#include "impl/particles_impl_sstp.ipp"
#include "impl/particles_impl_sstp_chem.ipp"
#include "impl/particles_impl_ante_adding_SD.ipp"
#include "impl/particles_impl_post_adding_SD.ipp"
#include "impl/particles_impl_src.ipp"
#include "impl/particles_impl_src_dry_distros_simple.ipp"
#include "impl/particles_impl_src_dry_distros_matching.ipp"
#include "impl/particles_impl_src_dry_distros.ipp"
#include "impl/particles_impl_src_dry_sizes.ipp"
#include "impl/particles_impl_rlx.ipp"
#include "impl/particles_impl_rlx_dry_distros.ipp"
#include "impl/particles_impl_update_incloud_time.ipp"
#include "impl/particles_impl_adjust_timesteps.ipp"

// MPI copy
#include "impl/particles_impl_pack.ipp"
#include "impl/particles_impl_unpack.ipp"
#include "impl/particles_impl_mpi_exchange.ipp"
