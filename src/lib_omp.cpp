#include "lib.hpp"

// workarounding Thrust bug #383 (Thanks to Jared for suggestion!)
#include <thrust/system/omp/execution_policy.h>

#include <thrust/system/omp/vector.h>
namespace thrust_device = ::thrust::omp;

#include "particles.tpp"
#include <omp.h>

namespace libcloudphxx
{ 
  namespace lgrngn
  {
    // checking if the above workaround actually did the job
    template <typename real_t, backend_t backend>
    void particles_t<real_t, backend>::impl::sanity_checks()
    {   
      if (omp_get_max_threads() == 1) return;

      thrust::omp::vector<int> v(100);

      struct 
      { 
        int operator()(int) const
        { 
          return omp_get_thread_num(); 
        } 
      } thread_id;

      thrust::transform(v.begin(), v.end(), v.begin(), thread_id);

      auto minmax = thrust::minmax_element(v.begin(), v.end());
      if (*minmax.first == *minmax.second)
        throw std::runtime_error("libcloudph++: OpenMP seems not to work properly!");
    }

    // instantiation 
    template class particles_t<float, OpenMP>;
    template class particles_t<double, OpenMP>;
  };
};
