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
        const arrinfo_t<real_t> courant_x
      ) { 
        this->init(th, rv, rhod, courant_x, arrinfo_t<real_t>(), arrinfo_t<real_t>()); 
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
      virtual void diag_sd_conc()                                   { assert(false); }
      virtual void diag_all()                                       { assert(false); }
      virtual void diag_rw_ge_rc()                                  { assert(false); }
      virtual void diag_RH_ge_Sc()                                  { assert(false); }
      virtual void diag_dry_rng(const real_t&, const real_t&)       { assert(false); }
      virtual void diag_wet_rng(const real_t&, const real_t&)       { assert(false); }
      virtual void diag_dry_mom(const int&)                         { assert(false); }
      virtual void diag_wet_mom(const int&)                         { assert(false); }
      virtual void diag_wet_mass_dens(const real_t&, const real_t&) { assert(false); }
      virtual void diag_chem(const enum chem_species_t&)            { assert(false); }
      virtual real_t *outbuf()                                      { assert(false); return NULL; }

      // storing a pointer to opts_init (e.g. for interrogatin about
      // dimensions in Python bindings)
      opts_init_t<real_t> *opts_init;

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
      void diag_wet_mass_dens(const real_t&, const real_t&);
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
      particles_t(const opts_init_t<real_t> &opts_init, const int &dev_id = -1, const int &n_cell_bfr = 0); // only opts_init specified by user

      // helper typedef
      typedef particles_proto_t<real_t> parent_t;
    };



    // specialization for the multi_GPU backend
    // has the init, stepping and diag functions
    // plus list of pointers to particles_t<CUDA> on each GPU
    // TODO: more elegant way?
    template <typename real_t>
    struct particles_t<real_t, multi_CUDA>: particles_proto_t<real_t>
    {
      // additional members
      std::vector<boost::shared_ptr<particles_t<real_t, CUDA > > > particles; // pointers to particles_t on each GPU
//      boost::ptr_vector<particles_t<real_t, CUDA> > particles; // pointer to particles_t on each GPU
#if defined(__NVCC__)
      cudaStream_t *streams; // cuda Stream on each device, used during P2P async memory copies
#endif
      opts_init_t<real_t> glob_opts_init; // global copy of opts_init (threads store their own in impl), 
      const int n_cell_tot;               // total number of cells
      std::vector<real_t> real_n_cell_tot; // vector of the size of the total number of cells to store output

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
      void diag_wet_mass_dens(const real_t&, const real_t&);
      real_t *outbuf();

      void diag_chem(const enum chem_species_t&);
      void diag_rw_ge_rc();
      void diag_RH_ge_Sc();
      void diag_all();

      // constructors
      particles_t(const opts_init_t<real_t> &opts_init, const int &dev_id = -1, const int &n_cell_bfr = 0); // only opts_init specified by user

      // helper typedef
      typedef particles_proto_t<real_t> parent_t;
    };
  };
};
