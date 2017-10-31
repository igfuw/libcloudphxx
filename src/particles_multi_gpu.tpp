// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief Thrust-based CPU/GPU particle-tracking logic for Lagrangian microphysics
  */

#include "detail/multi_gpu_utils.hpp"
#include "particles_multi_gpu_pimpl_ctor.ipp"
#include "particles_multi_gpu_diag.ipp"
#include "particles_multi_gpu_step.ipp"

