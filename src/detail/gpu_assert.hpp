#pragma once

#include <functional>
#include <thread>

#include <mutex>
#include <condition_variable>

#include <curand.h>

// macro to check for cuda errors, taken from 
// http://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
#define gpuErrchk(ans) { detail::gpuAssert((ans), __FILE__, __LINE__); }

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {   
      // https://stackoverflow.com/questions/40704147/built-in-function-to-get-curands-curandstatus-t-as-string
      #ifdef CURAND_H_
      // cuRAND API errors
      static const char *curandGetErrorString(curandStatus_t error)
      {
          switch (error)
          {
              case CURAND_STATUS_SUCCESS:
                  return "CURAND_STATUS_SUCCESS";
      
              case CURAND_STATUS_VERSION_MISMATCH:
                  return "CURAND_STATUS_VERSION_MISMATCH";
      
              case CURAND_STATUS_NOT_INITIALIZED:
                  return "CURAND_STATUS_NOT_INITIALIZED";
      
              case CURAND_STATUS_ALLOCATION_FAILED:
                  return "CURAND_STATUS_ALLOCATION_FAILED";
      
              case CURAND_STATUS_TYPE_ERROR:
                  return "CURAND_STATUS_TYPE_ERROR";
      
              case CURAND_STATUS_OUT_OF_RANGE:
                  return "CURAND_STATUS_OUT_OF_RANGE";
      
              case CURAND_STATUS_LENGTH_NOT_MULTIPLE:
                  return "CURAND_STATUS_LENGTH_NOT_MULTIPLE";
      
              case CURAND_STATUS_DOUBLE_PRECISION_REQUIRED:
                  return "CURAND_STATUS_DOUBLE_PRECISION_REQUIRED";
      
              case CURAND_STATUS_LAUNCH_FAILURE:
                  return "CURAND_STATUS_LAUNCH_FAILURE";
      
              case CURAND_STATUS_PREEXISTING_FAILURE:
                  return "CURAND_STATUS_PREEXISTING_FAILURE";
      
              case CURAND_STATUS_INITIALIZATION_FAILED:
                  return "CURAND_STATUS_INITIALIZATION_FAILED";
      
              case CURAND_STATUS_ARCH_MISMATCH:
                  return "CURAND_STATUS_ARCH_MISMATCH";
      
              case CURAND_STATUS_INTERNAL_ERROR:
                  return "CURAND_STATUS_INTERNAL_ERROR";
          }
      
          return "<unknown>";
      }
      #endif

      inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
      {   
         if (code != cudaSuccess) 
         {
            fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
            if (abort) exit(code);
         }
      }   

      inline void gpuAssert(curandStatus_t code, const char *file, int line, bool abort=true)
      {   
         if (code != CURAND_STATUS_SUCCESS) 
         {
            fprintf(stderr,"GPUassert (curand): %s %s %d\n", curandGetErrorString(code), file, line);
            if (abort) exit(code);
         }
      }   

      // max(1, n)
      inline int m1(int n) { return n == 0 ? 1 : n; }

      // run a function on a specific gpu
      inline void set_device_and_run(int id, std::function<void()> fun)
      {
        gpuErrchk(cudaSetDevice(id));
        fun();
      }

      // cxx_thread barrier
      // taken from libmpdata++, which in turn is based on boost barrier's code
      class barrier_t
      {
	std::mutex m_mutex;
	std::condition_variable m_cond;
	std::size_t m_generation, m_count;
        const std::size_t m_threshold;

	public:

	explicit barrier_t(const std::size_t count) : 
          m_count(count), 
          m_threshold(count),
          m_generation(0) 
        { }

	bool wait()
	{
          std::unique_lock<std::mutex> lock(m_mutex);
          unsigned int gen = m_generation;

          if (--m_count == 0)
          {
            m_generation++;
            m_count = m_threshold;
            m_cond.notify_all();
            return true;
          }

          while (gen == m_generation)
            m_cond.wait(lock);
          return false;
	}
      };
    }
  }
}
