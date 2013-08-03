#pragma once 

#include <cassert>
#include <memory>
#include <map>

#include "options.hpp"
#include "arrinfo.hpp"

namespace libcloudphxx
{
  namespace lgrngn
  {
    // to make inclusion of Thrust not neccesarry
    enum {cpp, omp, cuda}; 


    // to allow storing instances for multiple backends in one container/pointer
    template <typename real_t>
    struct particles_proto // TODO: rename to any?
    {
      // dataZero + strides

      // 3D version
      virtual void init(
        const arrinfo_t<real_t> rhod_th,
        const arrinfo_t<real_t> rhod_rv,
        const arrinfo_t<real_t> rhod, 
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_y,
        const arrinfo_t<real_t> courant_z
      ) { assert(false); }  

      // 2D version
      void init(
        const arrinfo_t<real_t> rhod_th,
        const arrinfo_t<real_t> rhod_rv,
        const arrinfo_t<real_t> rhod,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_z
      ) { this->init(rhod_th, rhod_rv, rhod, courant_x, arrinfo_t<real_t>(), courant_z); }  
 
      // 1D version
      void init(
        const arrinfo_t<real_t> rhod_th,
        const arrinfo_t<real_t> rhod_rv,
        const arrinfo_t<real_t> rhod,
        const arrinfo_t<real_t> courant_z
      ) { this->init(rhod_th, rhod_rv, rhod, arrinfo_t<real_t>(), arrinfo_t<real_t>(), courant_z); }  

      // 0D version
      void init(
        const arrinfo_t<real_t> rhod_th,
        const arrinfo_t<real_t> rhod_rv,
        const arrinfo_t<real_t> rhod
      ) { this->init(rhod_th, rhod_th, rhod, arrinfo_t<real_t>(), arrinfo_t<real_t>(), arrinfo_t<real_t>()); }  




      // 3D variable density version
      virtual void step(
        arrinfo_t<real_t> rhod_th,
        arrinfo_t<real_t> rhod_rv,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_y,
        const arrinfo_t<real_t> courant_z,
        const arrinfo_t<real_t> rhod
      ) { assert(false); }  

      // 3D constant density version
      void step(
        arrinfo_t<real_t> rhod_th,
        arrinfo_t<real_t> rhod_rv,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_y,
        const arrinfo_t<real_t> courant_z
      ) { this->step(rhod_th, rhod_rv, courant_x, courant_y, courant_z, arrinfo_t<real_t>()); }  

      // 2D constant density version
      void step(
        arrinfo_t<real_t> rhod_th,
        arrinfo_t<real_t> rhod_rv,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_z
      ) { this->step(rhod_th, rhod_rv, courant_x, arrinfo_t<real_t>(), courant_z, arrinfo_t<real_t>()); }  

      // 1D constant density version
      void step(
        arrinfo_t<real_t> rhod_th,
        arrinfo_t<real_t> rhod_rv,
        const arrinfo_t<real_t> courant_z
      ) { this->step(rhod_th, rhod_rv, arrinfo_t<real_t>(), arrinfo_t<real_t>(), courant_z, arrinfo_t<real_t>()); }  

      // 0D constant density version
      void step(
        arrinfo_t<real_t> rhod_th,
        arrinfo_t<real_t> rhod_rv
      ) { this->step(rhod_th, rhod_rv, arrinfo_t<real_t>(), arrinfo_t<real_t>(), arrinfo_t<real_t>(), arrinfo_t<real_t>()); }  

      virtual void diag(
      ) { assert(false); }

      virtual real_t *outbuf() 
      { 
        assert(false);
        return NULL;
      }
    };  

    // prototype of what's implemented in the .tpp file
    template <typename real_t, int thrust_device_system>
    class particles : public particles_proto<real_t>
    {
      typedef particles_proto<real_t> parent_t;

      // pimpl stuff
      struct impl;
      std::auto_ptr<impl> pimpl;
  
      // the public API
      public:  
      particles(const opts_t<real_t>); // ctor

      // init separated from the ctor as not all data might be available
      void init(
        const arrinfo_t<real_t> rhod_th,
        const arrinfo_t<real_t> rhod_rv,
        const arrinfo_t<real_t> rhod,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_y, 
        const arrinfo_t<real_t> courant_z
      );

      void step(
        const arrinfo_t<real_t> rhod_th,
        const arrinfo_t<real_t> rhod_rv,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_y,
        const arrinfo_t<real_t> courant_z,
        const arrinfo_t<real_t> rhod 
      );

      void diag();

      real_t *outbuf();
    };

    // to be explicitely instantiated
    template <typename real_t>
    struct factory
    {
      static particles_proto<real_t> *make(const int backend, const opts_t<real_t> &);
    };
  };
};
