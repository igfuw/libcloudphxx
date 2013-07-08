// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief Thrust-based CPU/GPU particle-tracking logic for Lagrangian microphysics
  */

#pragma once

#include <iostream>

#include "particles.hpp"

#include <thrust/device_vector.h>

namespace libcloudphxx
{
  namespace common
  {
    namespace prtcls
    {

      template <typename real_t, int device>
      void particles<real_t, device>::func()
      {
#if (THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_OMP) 
	assert(device == omp);
#elif (THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA) 
	assert(device == cuda);
#elif (THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CPP) 
	assert(device == cpp);
#else
	assert(false);
#endif

	std::cerr << "CUDA/OpenMP/CPP: " << device << std::endl;
	thrust::device_vector<real_t> vec(1024*1024);
      }
    };
  };
};
