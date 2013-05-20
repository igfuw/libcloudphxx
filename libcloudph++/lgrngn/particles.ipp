/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  *
  * @brief Thrust-based CPU/GPU particle-tracking logic for Lagrangian microphysics
  * @details this file is to included in a .cpp file! For example:
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

#include <iostream>

#include <libcloudph++/lgrngn/particles.hpp>

#include <thrust/device_vector.h>

template <typename real_t, int thrust_device_system>
void particles<real_t, thrust_device_system>::func()
{
#if (THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_OMP) 
  assert(device_system_macro == openmp);
#elif (THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA) 
  assert(device_system_macro == cuda);
#elif (THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CPP) 
  assert(device_system_macro == cpp);
#endif

  std::cerr << "CUDA/OpenMP/CPP: " << thrust_device_system << std::endl;
  thrust::device_vector<real_t> vec(1024*1024);
}

template class particles<float, device_system_macro>;
