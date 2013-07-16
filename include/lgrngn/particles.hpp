#pragma once 

#include <cassert>
#include <memory>

#include "../common/unary_function.hpp"

namespace libcloudphxx
{
  namespace lgrngn
  {
    using libcloudphxx::common::unary_function;

    // to make inclusion of Thrust not neccesarry
    enum {cpp, omp, cuda};

    // to allow storing instances for multiple backends in one container/pointer
    template <typename real_t>
    class particles_proto // TODO: rename to any?
    {
      public: 
      virtual void init(unary_function<real_t> *pdf) { assert(false); }  

      // TODO: -> pimpl? / private?
      virtual void init_dry(unary_function<real_t> *pdf) { assert(false); }  
      virtual void init_xyz() { assert(false); }  
      virtual void init_Tpr() { assert(false); }  
      virtual void init_wet() { assert(false); }  
    };  

    template <typename real_t>
    struct opts_t 
    {   
      int nx, ny, nz;
      real_t sd_conc_mean; 
      real_t dx, dy, dz;

      // defaults
      opts_t() 
      { 
        nx = ny = nz = 0; 
        dx = dy = dz = 1;
        sd_conc_mean = -1;
      }
    };  

    // prototype of what's implemented in the .tpp file
    template <typename real_t, int thrust_device_system>
    class particles : public particles_proto<real_t>
    {
      // pimpl stuff
      struct impl;
      std::auto_ptr<impl> pimpl;
  
      // the public API
      public:  
      particles(const opts_t<real_t>); // ctor

      void init(unary_function<real_t> *pdf);
      // TODO: -> pimpl? / private?
      void init_dry(unary_function<real_t> *pdf);
      void init_xyz();
      void init_Tpr();
      void init_wet();
    };

    // to be explicitely instantiated
    template <typename real_t>
    particles_proto<real_t> *factory(const int backend, const real_t sd_conc_mean); // 0D version

    template <typename real_t>
    particles_proto<real_t> *factory(const int backend, const real_t sd_conc_mean, 
      const int nz, const real_t dz); // 1D version

    template <typename real_t>
    particles_proto<real_t> *factory(const int backend, const real_t sd_conc_mean,
      const int, const int, const real_t, const real_t); // 2D version

    template <typename real_t>
    particles_proto<real_t> *factory(const int backend, const real_t sd_conc_mean,
      const int, const int, const int, const real_t, const real_t, const real_t); // 3D version
  };
};
