// including it first not to require pthread option to nvcc
//#include <blitz/array.h>

// TODO: setting OMP in gcc-compiled binaries seems needed
//       both to get OpenMP working on a CUDA system (rintime)
//       and to get it compiled on a non-CUDA system
//       setting it here however makes nvcc complain
//       but as of now there's a conficlict of definition of default::... TODO!
//#define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_OMP
#include <thrust/system/cuda/vector.h>
namespace thrust_device = ::thrust::cuda;

#include "particles.tpp"

namespace libcloudphxx
{ 
  namespace lgrngn
  {
    // instantiation 
    template class particles<float, cuda>;
//    template class particles<double, cuda>;
  };
};
