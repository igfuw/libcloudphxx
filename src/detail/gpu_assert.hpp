#pragma once

// macro to check for cuda errors, taken from 
// move it to utils...
// http://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
#define gpuErrchk(ans) { detail::gpuAssert((ans), __FILE__, __LINE__); }

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {   
      void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
      {   
         if (code != cudaSuccess) 
         {
            fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
            if (abort) exit(code);
         }
      }   

      // max(1, n)
      int m1(int n) { return n == 0 ? 1 : n; }
    }
  }
}
