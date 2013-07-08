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
      typedef thrust::device_vector<int>::size_type thrust_size_t;

      // pimpl stuff
      template <typename real_t, int device>
      struct particles<real_t, device>::detail
      { 
        const int n_dims;
        thrust_size_t n_part;
        thrust::device_vector<real_t> 
          rd3, // dry radii cubed 
          xi, 
          x, 
          y, 
          z;

        //
        void init()
        {
	  // sanity checks (Thrust preprocessor macro names vs. local enum names)
#if (THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_OMP) 
	  assert(device == omp);
#elif (THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CUDA) 
	  assert(device == cuda);
#elif (THRUST_DEVICE_SYSTEM == THRUST_DEVICE_SYSTEM_CPP) 
	  assert(device == cpp);
#else
	  assert(false);
#endif

          // 
          rd3.resize(n_part);
        }

        // 0D ctor 
        detail(real_t sd_conc_mean)
          : n_dims(0), n_part(sd_conc_mean)
        {
std::cerr << "aqq" << std::endl;
          init();
        }
      };

      // 0D ctor
      template <typename real_t, int device>
      particles<real_t, device>::particles(real_t sd_conc_mean)
        : pimpl(new detail(sd_conc_mean)) 
      {}
    };
  };
};
