// including it first not to require pthread option to nvcc
//#include <blitz/array.h>

#define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_CUDA
#define libcloudphxx_particles_device cuda
#define libcloudphxx_particles_real_t float
#include "particles.ipp"
