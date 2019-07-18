#pragma once
namespace
{
  #if !defined(__NVCC__)
    using std::isnan;
    using std::isinf;
  #endif
  
  // functor for testing nan/inf values
  struct isnaninf { 
    template<class real_t>
    BOOST_GPU_ENABLED
    bool operator()(const real_t a) const {
      return isnan(a) || isinf(a);
    }
  };
};

#ifdef NDEBUG
#define nancheck(arr, name) ((void)0)
#else
#define nancheck(arr, name) {\
  int nan_count = thrust::transform_reduce(arr.begin(), arr.end(), isnaninf(), 0, thrust::plus<bool>());\
  if(nan_count>0){\
    std::cout << nan_count << " nan/inf numbers detected in: " << name << std::endl;\
    debug::print(arr);\
    assert(0);}}
#endif

#ifdef NDEBUG
#define nancheck_range(begin, end, name) ((void)0)
#else
#define nancheck_range(begin, end, name) {\
  int nan_count = thrust::transform_reduce(begin, end, isnaninf(), 0, thrust::plus<bool>());\
  if(nan_count>0){\
    std::cout << nan_count << " nan/inf numbers detected in: " << name << std::endl;\
    debug::print(begin, end);\
    assert(0);}}
#endif
