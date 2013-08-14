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
    enum { cpp, omp, cuda }; 

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
      virtual void step_sync(
        arrinfo_t<real_t> rhod_th,
        arrinfo_t<real_t> rhod_rv,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_y,
        const arrinfo_t<real_t> courant_z,
        const arrinfo_t<real_t> rhod
      ) { assert(false); }  

      virtual void step_async() { assert(false); }  

      // 3D constant density version
      void step_sync(
        arrinfo_t<real_t> rhod_th,
        arrinfo_t<real_t> rhod_rv,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_y,
        const arrinfo_t<real_t> courant_z
      ) { this->step_sync(rhod_th, rhod_rv, courant_x, courant_y, courant_z, arrinfo_t<real_t>()); }  

      // 2D constant density version
      void step_sync(
        arrinfo_t<real_t> rhod_th,
        arrinfo_t<real_t> rhod_rv,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_z
      ) { this->step_sync(rhod_th, rhod_rv, courant_x, arrinfo_t<real_t>(), courant_z, arrinfo_t<real_t>()); }  

      // 1D constant density version
      void step_sync(
        arrinfo_t<real_t> rhod_th,
        arrinfo_t<real_t> rhod_rv,
        const arrinfo_t<real_t> courant_z
      ) { this->step_sync(rhod_th, rhod_rv, arrinfo_t<real_t>(), arrinfo_t<real_t>(), courant_z, arrinfo_t<real_t>()); }  

      // 0D constant density version
      void step_sync(
        arrinfo_t<real_t> rhod_th,
        arrinfo_t<real_t> rhod_rv
      ) { this->step_sync(rhod_th, rhod_rv, arrinfo_t<real_t>(), arrinfo_t<real_t>(), arrinfo_t<real_t>(), arrinfo_t<real_t>()); }  

      // method for accessing super-droplet statistics
      virtual void diag_sd_conc()                             { assert(false); }
      virtual void diag_dry_rng(const real_t&, const real_t&) { assert(false); }
      virtual void diag_wet_rng(const real_t&, const real_t&) { assert(false); }
      virtual void diag_dry_mom(const int&)                   { assert(false); }
      virtual void diag_wet_mom(const int&)                   { assert(false); }
      virtual real_t *outbuf()                                { assert(false); return NULL; }
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
      particles(const opts_t<real_t> &); // ctor

      // init separated from the ctor as not all data might be available
      void init(
        const arrinfo_t<real_t> rhod_th,
        const arrinfo_t<real_t> rhod_rv,
        const arrinfo_t<real_t> rhod,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_y, 
        const arrinfo_t<real_t> courant_z
      );

      void step_sync(
        arrinfo_t<real_t> rhod_th,
        arrinfo_t<real_t> rhod_rv,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_y,
        const arrinfo_t<real_t> courant_z,
        const arrinfo_t<real_t> rhod 
      );

      void step_async();

      void diag_sd_conc();
      void diag_dry_rng(const real_t&, const real_t&);
      void diag_wet_rng(const real_t&, const real_t&);
      void diag_dry_mom(const int&);
      void diag_wet_mom(const int&);
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
