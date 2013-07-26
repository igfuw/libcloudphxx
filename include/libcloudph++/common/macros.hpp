#pragma once

#if __cplusplus < 201103L
//#  warning using workarounds instead of C++11 construct
#  define constexpr const
#endif

#define libcloudphxx_const(type, name, value, unit) \
  template <typename real_t> \
  BOOST_GPU_ENABLED \
  static constexpr quantity<type, real_t> name() { return real_t(value) * unit; }

#define libcloudphxx_const_derived(type, name, value) \
  template <typename real_t> \
  BOOST_GPU_ENABLED \
  static constexpr quantity<type, real_t> name() { return value; }
