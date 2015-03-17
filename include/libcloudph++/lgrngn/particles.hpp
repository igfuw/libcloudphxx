#pragma once 

#include <libcloudph++/lgrngn/extincl.hpp>

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
        const arrinfo_t<real_t> th,
        const arrinfo_t<real_t> rv,
        const arrinfo_t<real_t> rhod, 
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_y,
        const arrinfo_t<real_t> courant_z
      ) { 
        assert(false); 
      }  

      // 2D version
      void init(
        const arrinfo_t<real_t> th,
        const arrinfo_t<real_t> rv,
        const arrinfo_t<real_t> rhod,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_z
      ) { 
        this->init(th, rv, rhod, courant_x, arrinfo_t<real_t>(), courant_z); 
      }  
 
      // 1D version
      void init(
        const arrinfo_t<real_t> th,
        const arrinfo_t<real_t> rv,
        const arrinfo_t<real_t> rhod,
        const arrinfo_t<real_t> courant_z
      ) { 
        this->init(th, rv, rhod, arrinfo_t<real_t>(), arrinfo_t<real_t>(), courant_z); 
      }  

      // 0D version
      void init(
        const arrinfo_t<real_t> th,
        const arrinfo_t<real_t> rv,
        const arrinfo_t<real_t> rhod
      ) { 
        this->init(th, rv, rhod, arrinfo_t<real_t>(), arrinfo_t<real_t>(), arrinfo_t<real_t>()); 
      }  

      // 3D variable density version
      virtual void step_sync(
        const opts_t<real_t> &,
        arrinfo_t<real_t> th,
        arrinfo_t<real_t> rv,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_y,
        const arrinfo_t<real_t> courant_z,
        const arrinfo_t<real_t> rhod
      ) { 
        assert(false); 
      }  

      // returns accumulated rainfall
      virtual real_t step_async(
        const opts_t<real_t> &
      ) { 
        assert(false); 
        return 0;
      }  

      // 3D constant density version
      void step_sync(
        const opts_t<real_t> &opts,
        arrinfo_t<real_t> th,
        arrinfo_t<real_t> rv,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_y,
        const arrinfo_t<real_t> courant_z
      ) { 
        this->step_sync(opts, th, rv, courant_x, courant_y, courant_z, arrinfo_t<real_t>()); 
      }  

      // kinematic version
      void step_sync(
        const opts_t<real_t> &opts,
        arrinfo_t<real_t> th,
        arrinfo_t<real_t> rv
      ) { 
        this->step_sync(
          opts,
          th, 
          rv, 
          arrinfo_t<real_t>(), 
          arrinfo_t<real_t>(), 
          arrinfo_t<real_t>(), 
          arrinfo_t<real_t>()
        ); 
      }  


      // 2D constant density version
      void step_sync(
        const opts_t<real_t> &opts,
        arrinfo_t<real_t> th,
        arrinfo_t<real_t> rv,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_z
      ) { 
        this->step_sync(
          opts,
          th, 
          rv, 
          courant_x, 
          arrinfo_t<real_t>(), 
          courant_z, 
          arrinfo_t<real_t>()
        ); 
      }  

      // 0D version
      void step_sync(
        const opts_t<real_t> &opts,
        arrinfo_t<real_t> th,
        arrinfo_t<real_t> rv,
        const arrinfo_t<real_t> rhod
      ) { 
        this->step_sync(
          opts,
          th, 
          rv, 
          arrinfo_t<real_t>(), 
          arrinfo_t<real_t>(), 
          arrinfo_t<real_t>(), 
          rhod
        ); 
      }  

      // method for accessing super-droplet statistics
      virtual void diag_sd_conc()                             { assert(false); }
      virtual void diag_all()                                 { assert(false); }
      virtual void diag_rw_ge_rc()                            { assert(false); }
      virtual void diag_RH_ge_Sc()                            { assert(false); }
      virtual void diag_dry_rng(const real_t&, const real_t&) { assert(false); }
      virtual void diag_wet_rng(const real_t&, const real_t&) { assert(false); }
      virtual void diag_dry_mom(const int&)                   { assert(false); }
      virtual void diag_wet_mom(const int&)                   { assert(false); }
      virtual void diag_chem(const enum chem_species_t&)      { assert(false); }
      virtual real_t *outbuf()                                { assert(false); return NULL; }

      // storing a pointer to opts_init (e.g. for interrogatin about
      // dimensions in Python bindings)
      const opts_init_t<real_t> *opts_init;
    };  

    // prototype of what's implemented in the .tpp file
//<listing>
    template <typename real_t, backend_t backend>
    struct particles_t: particles_proto_t<real_t>
    {
      // initialisation 
      void init(
        const arrinfo_t<real_t> th,
        const arrinfo_t<real_t> rv,
        const arrinfo_t<real_t> rhod,
        const arrinfo_t<real_t> courant_1,
        const arrinfo_t<real_t> courant_2, 
        const arrinfo_t<real_t> courant_3
      );

      // time-stepping methods
      void step_sync(
        const opts_t<real_t> &,
        arrinfo_t<real_t> th,
        arrinfo_t<real_t> rv,
        const arrinfo_t<real_t> courant_1,
        const arrinfo_t<real_t> courant_2,
        const arrinfo_t<real_t> courant_3,
        const arrinfo_t<real_t> rhod 
      );
      real_t step_async(
        const opts_t<real_t> &
      );

      // diagnostic methods
      void diag_sd_conc();
      void diag_dry_rng(
        const real_t &r_mi, const real_t &r_mx
      );
      void diag_wet_rng(
        const real_t &r_mi, const real_t &r_mx
      );
      void diag_dry_mom(const int &k);
      void diag_wet_mom(const int &k);
      real_t *outbuf();

      // ...
//</listing>

      void diag_chem(const enum chem_species_t&);
      void diag_rw_ge_rc();
      void diag_RH_ge_Sc();
      void diag_all();

      struct impl;
      std::auto_ptr<impl> pimpl;

      // constructor
      particles_t(const opts_init_t<real_t> &);

      // helper typedef
      typedef particles_proto_t<real_t> parent_t;
    };
  };
};
