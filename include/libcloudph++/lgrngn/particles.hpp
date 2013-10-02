#pragma once 

#include <cassert>
#include <memory>
#include <map>

#include "opts.hpp"
#include "opts_init.hpp"
#include "arrinfo.hpp"
#include "backend.hpp"

namespace libcloudphxx
{
  namespace lgrngn
  {
    // to allow storing instances for multiple backends in one container/pointer
    template <typename real_t>
    struct particles_proto_t 
    {
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
        const opts_t<real_t> &,
        arrinfo_t<real_t> rhod_th,
        arrinfo_t<real_t> rhod_rv,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_y,
        const arrinfo_t<real_t> courant_z,
        const arrinfo_t<real_t> rhod
      ) { assert(false); }  

      virtual void step_async(
        const opts_t<real_t> &
      ) { assert(false); }  

      // 3D constant density version
      void step_sync(
        const opts_t<real_t> &opts,
        arrinfo_t<real_t> rhod_th,
        arrinfo_t<real_t> rhod_rv,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_y,
        const arrinfo_t<real_t> courant_z
      ) { this->step_sync(opts, rhod_th, rhod_rv, courant_x, courant_y, courant_z, arrinfo_t<real_t>()); }  

      // 2D constant density version
      void step_sync(
        const opts_t<real_t> &opts,
        arrinfo_t<real_t> rhod_th,
        arrinfo_t<real_t> rhod_rv,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_z
      ) { 
        this->step_sync(
          opts,
          rhod_th, 
          rhod_rv, 
          courant_x, 
          arrinfo_t<real_t>(), 
          courant_z, 
          arrinfo_t<real_t>()
        ); 
      }  

      // 1D constant density version
      void step_sync(
        const opts_t<real_t> &opts,
        arrinfo_t<real_t> rhod_th,
        arrinfo_t<real_t> rhod_rv,
        const arrinfo_t<real_t> courant_z
      ) { 
        this->step_sync(
          opts,
          rhod_th, 
          rhod_rv, 
          arrinfo_t<real_t>(),
          arrinfo_t<real_t>(), 
          courant_z, 
          arrinfo_t<real_t>()); 
      }  

      // 0D constant density version
      void step_sync(
        const opts_t<real_t> &opts,
        arrinfo_t<real_t> rhod_th,
        arrinfo_t<real_t> rhod_rv
      ) { 
        this->step_sync(
          opts,
          rhod_th, 
          rhod_rv, 
          arrinfo_t<real_t>(), 
          arrinfo_t<real_t>(), 
          arrinfo_t<real_t>(), 
          arrinfo_t<real_t>()
        ); 
      }  

      // method for accessing super-droplet statistics
      virtual void diag_sd_conc()                             { assert(false); }
      virtual void diag_dry_rng(const real_t&, const real_t&) { assert(false); }
      virtual void diag_wet_rng(const real_t&, const real_t&) { assert(false); }
      virtual void diag_dry_mom(const int&)                   { assert(false); }
      virtual void diag_wet_mom(const int&)                   { assert(false); }
      virtual real_t *outbuf()                                { assert(false); return NULL; }
    };  

    // prototype of what's implemented in the .tpp file
//<listing>
    template <typename real_t, backend_t backend>
    struct particles_t: particles_proto_t<real_t>
    {
      // initialisation 
      void init(
        const arrinfo_t<real_t> rho_e,
        const arrinfo_t<real_t> rho_v,
        const arrinfo_t<real_t> rho_d,
        const arrinfo_t<real_t> courant_1,
        const arrinfo_t<real_t> courant_2, 
        const arrinfo_t<real_t> courant_3
      );

      // time-stepping methods
      void step_sync(
        const opts_t<real_t> &,
        arrinfo_t<real_t> rho_e,
        arrinfo_t<real_t> rho_v,
        const arrinfo_t<real_t> courant_1,
        const arrinfo_t<real_t> courant_2,
        const arrinfo_t<real_t> courant_3,
        const arrinfo_t<real_t> rho_d 
      );
      void step_async(
        const opts_t<real_t> &
      );

      // diagnostic methods
      void diag_sd_conc();
      void diag_dry_rng(
        const real_t &min, const real_t &max
      );
      void diag_wet_rng(
        const real_t &min, const real_t &max
      );
      void diag_dry_mom(const int &num);
      void diag_wet_mom(const int &num);
      real_t *outbuf();

      // ...
//</listing>
      struct impl;
      std::auto_ptr<impl> pimpl;

      // constructor
      particles_t(const opts_init_t<real_t> &);

      // helper typedef
      typedef particles_proto_t<real_t> parent_t;
    };
  };
};
