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
      // initialisation 
      virtual void init(
        const arrinfo_t<real_t> th,
        const arrinfo_t<real_t> rv,
        const arrinfo_t<real_t> rhod,
        const arrinfo_t<real_t> courant_x = arrinfo_t<real_t>(),
        const arrinfo_t<real_t> courant_y = arrinfo_t<real_t>(), 
        const arrinfo_t<real_t> courant_z = arrinfo_t<real_t>(),
        const std::map<enum chem_species_t, const arrinfo_t<real_t> > ambient_chem = std::map<enum chem_species_t, const arrinfo_t<real_t> >()
      ) { 
        assert(false);
      }  
 
      // stuff that requires Eulerian component to wait
      virtual void step_sync(
        const opts_t<real_t> &,
        arrinfo_t<real_t> th,
        arrinfo_t<real_t> rv,
        const arrinfo_t<real_t> rhod      = arrinfo_t<real_t>(),
        const arrinfo_t<real_t> courant_x = arrinfo_t<real_t>(),
        const arrinfo_t<real_t> courant_y = arrinfo_t<real_t>(),
        const arrinfo_t<real_t> courant_z = arrinfo_t<real_t>(),
        std::map<enum chem_species_t, arrinfo_t<real_t> > ambient_chem = std::map<enum chem_species_t, arrinfo_t<real_t> >()
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

      // add some cloud water
      virtual void step_rc_adjust(
        const arrinfo_t<real_t> rc_adj
      ) { 
        assert(false); 
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
      virtual void diag_precip_rate()                               { assert(false); }
      virtual real_t *outbuf()                                      { assert(false); return NULL; }

      // storing a pointer to opts_init (e.g. for interrogatin about
      // dimensions in Python bindings)
      opts_init_t<real_t> *opts_init;

    };  

    // prototype of what's implemented in the .tpp file
    template <typename real_t, backend_t backend>
    struct particles_t: particles_proto_t<real_t>
    {
      // initialisation 
      void init(
        const arrinfo_t<real_t> th,
        const arrinfo_t<real_t> rv,
        const arrinfo_t<real_t> rhod,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_y, 
        const arrinfo_t<real_t> courant_z,
        const std::map<enum chem_species_t, const arrinfo_t<real_t> > ambient_chem
      );

      // time-stepping methods
      void step_sync(
        const opts_t<real_t> &,
        arrinfo_t<real_t> th,
        arrinfo_t<real_t> rv,
        const arrinfo_t<real_t> rhod,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_y,
        const arrinfo_t<real_t> courant_z,
        std::map<enum chem_species_t, arrinfo_t<real_t> > ambient_chem
      );

      real_t step_async(
        const opts_t<real_t> &
      );

      void step_rc_adjust(
        const arrinfo_t<real_t> 
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

      void diag_chem(const enum chem_species_t&);
      void diag_rw_ge_rc();
      void diag_RH_ge_Sc();
      void diag_all();
      void diag_precip_rate();

      struct impl;
      std::auto_ptr<impl> pimpl;

      // constructor
      particles_t(const opts_init_t<real_t> &opts_init, const int &n_cell_bfr = 0); // only opts_init specified by user

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
      boost::ptr_vector<particles_t<real_t, CUDA> > particles; // pointer to particles_t on each GPU
      opts_init_t<real_t> glob_opts_init; // global copy of opts_init (threads store their own in impl), 
      const int n_cell_tot;               // total number of cells
      std::vector<real_t> real_n_cell_tot; // vector of the size of the total number of cells to store output

      // initialisation 
      void init(
        const arrinfo_t<real_t> th,
        const arrinfo_t<real_t> rv,
        const arrinfo_t<real_t> rhod,
        const arrinfo_t<real_t> courant_x = arrinfo_t<real_t>(),
        const arrinfo_t<real_t> courant_y = arrinfo_t<real_t>(), 
        const arrinfo_t<real_t> courant_z = arrinfo_t<real_t>(),
        const std::map<enum chem_species_t, const arrinfo_t<real_t> > ambient_chem = std::map<enum chem_species_t, const arrinfo_t<real_t> >()
      );

      // time-stepping methods
      void step_sync(
        const opts_t<real_t> &,
        arrinfo_t<real_t> th,
        arrinfo_t<real_t> rv,
        const arrinfo_t<real_t> rhod      = arrinfo_t<real_t>(),
        const arrinfo_t<real_t> courant_x = arrinfo_t<real_t>(),
        const arrinfo_t<real_t> courant_y = arrinfo_t<real_t>(),
        const arrinfo_t<real_t> courant_z = arrinfo_t<real_t>(),
        std::map<enum chem_species_t, arrinfo_t<real_t> > ambient_chem = std::map<enum chem_species_t, arrinfo_t<real_t> >()
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
      void diag_precip_rate();

      // constructors
      particles_t(const opts_init_t<real_t> &opts_init, const int &n_x_bfr = 0); // only opts_init specified by user

      // helper typedef
      typedef particles_proto_t<real_t> parent_t;
    };
  };
};
