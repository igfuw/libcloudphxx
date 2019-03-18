#pragma once

// functor for testing not-positive values
struct isneg { 
  template<class real_t>
  BOOST_GPU_ENABLED
  bool operator()(const real_t a) const {
    return a <= real_t(0);
  }
};

#define _LIBCLOUDPHXX_SMALL_POSITIVE_VALUE 1e-10

// make not-positive values slightly positive...
#ifdef NDEBUG
#define negtozero(arr, name) {\
  thrust::replace_if(\
    arr.begin(), arr.end(), \
    isneg(), \
    real_t(_LIBCLOUDPHXX_SMALL_POSITIVE_VALUE) \
  );}
#else
#define negtozero(arr, name) {\
  if(*thrust::min_element(arr.begin(), arr.end()) <= 0.)\
  {\
    printf("A non-positive number detected in: " name "\n");\
    debug::print(arr);\
    printf("Cheating in libcloudphxx: replacing non-positive values with small positive ones\n");\
    thrust::replace_if(\
      arr.begin(), arr.end(), \
      isneg(), \
      real_t(_LIBCLOUDPHXX_SMALL_POSITIVE_VALUE) \
    );}}
#endif
