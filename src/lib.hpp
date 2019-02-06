#pragma once

// workarounding Thrust bug #382 (this is just a default setting not used hereinafter)
#if !defined(__NVCC__)
  #define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_CPP
#endif
