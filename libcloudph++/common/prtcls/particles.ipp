/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  *
  * @brief Thrust-based CPU/GPU particle-tracking logic for Lagrangian microphysics
  * @details this file is to be included in a .cpp file! For example:
  * 
  * a cpp file with a CUDA instance (to be conpiled with nvcc):
  * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  * #define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_CUDA
  * #define device_system_macro cuda
  * #include <libcloudph++/lgrngn/particles.ipp>
  * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  * 
  * and another cpp file making use of this instance (can be compiled with g++ or clang):
  * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  * ...
  * #include <libcloudph++/lgrngn/particles.hpp>
  * ...
  * p = new particles<float, cuda>();
  * ...
  * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  */

#include "particles.tpp"

// instantiation 
template class particles<
  libcloudphxx_particles_real_t, // float, double, ... 
  libcloudphxx_particles_device  // openmp, cuda, cpp
>;
