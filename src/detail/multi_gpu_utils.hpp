#pragma once

#include <functional>
#include <thread>

// macro to check for cuda errors, taken from 
// http://stackoverflow.com/questions/14038589/what-is-the-canonical-way-to-check-for-errors-using-the-cuda-runtime-api
#define gpuErrchk(ans) { detail::gpuAssert((ans), __FILE__, __LINE__); }

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {   
      template<class real_t>
      int get_dev_nx(const opts_init_t<real_t> &opts_init, const int &dev_no)
      {
        if(dev_no < opts_init.dev_count-1)
          return opts_init.nx / opts_init.dev_count + .5;
        else
          return opts_init.nx - dev_no * int(opts_init.nx / opts_init.dev_count + .5); 
      } 

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

      // run a function on a specific gpu
      void set_device_and_run(int id, std::function<void()> fun)
      {
        gpuErrchk(cudaSetDevice(id));
        fun();
      }

/*
      // run a function concurently on gpus, TODO: some forwarding/moving/passing by reference?
      void mcuda_run(int dev_count, std::function<void()> fun)
      {
        std::vector<std::thread> threads;
        for (int i = 0; i < dev_count; ++i)
        {
          threads.emplace_back(
            detail::set_device_and_run, i, fun      
          );
        }
        for (auto &th : threads) th.join();
      }
*/
    }
  }
}
