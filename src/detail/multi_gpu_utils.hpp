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
