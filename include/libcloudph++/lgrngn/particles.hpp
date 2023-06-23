#pragma once 

#include "extincl.hpp"

#include "opts.hpp"
#include "../common/output.hpp"
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
        const arrinfo_t<real_t> p = arrinfo_t<real_t>(),
        const arrinfo_t<real_t> courant_x = arrinfo_t<real_t>(),
        const arrinfo_t<real_t> courant_y = arrinfo_t<real_t>(), 
        const arrinfo_t<real_t> courant_z = arrinfo_t<real_t>(),
        const std::map<enum common::chem::chem_species_t, const arrinfo_t<real_t> > ambient_chem = std::map<enum common::chem::chem_species_t, const arrinfo_t<real_t> >()
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
        const arrinfo_t<real_t> diss_rate = arrinfo_t<real_t>(), // TKE dissipation rate (epsilon)
        std::map<enum common::chem::chem_species_t, arrinfo_t<real_t> > ambient_chem = std::map<enum common::chem::chem_species_t, arrinfo_t<real_t> >()
      ) { 
        assert(false); 
      }  

      virtual void sync_in(
        arrinfo_t<real_t> th,
        arrinfo_t<real_t> rv,
        const arrinfo_t<real_t> rhod      = arrinfo_t<real_t>(),
        const arrinfo_t<real_t> courant_x = arrinfo_t<real_t>(),
        const arrinfo_t<real_t> courant_y = arrinfo_t<real_t>(),
        const arrinfo_t<real_t> courant_z = arrinfo_t<real_t>(),
        const arrinfo_t<real_t> diss_rate = arrinfo_t<real_t>(),
        std::map<enum common::chem::chem_species_t, arrinfo_t<real_t> > ambient_chem = std::map<enum common::chem::chem_species_t, arrinfo_t<real_t> >()
      ) { 
        assert(false); 
      }  

      virtual void step_cond(
        const opts_t<real_t> &,
        arrinfo_t<real_t> th,
        arrinfo_t<real_t> rv,
        std::map<enum common::chem::chem_species_t, arrinfo_t<real_t> > ambient_chem = std::map<enum common::chem::chem_species_t, arrinfo_t<real_t> >()
      ) { 
        assert(false); 
      }  

      // returns accumulated rainfall
      virtual void step_async(
        const opts_t<real_t> &
      ) { 
        assert(false); 
      }  

      // method for accessing super-droplet statistics
      virtual void diag_sd_conc()                                               { assert(false); }
      virtual void diag_pressure()                                              { assert(false); }
      virtual void diag_temperature()                                           { assert(false); }
      virtual void diag_RH()                                                    { assert(false); }
      virtual void diag_all()                                                   { assert(false); }
      virtual void diag_rw_ge_rc()                                              { assert(false); }
      virtual void diag_RH_ge_Sc()                                              { assert(false); }
      virtual void diag_dry_rng(const real_t&, const real_t&)                   { assert(false); }
      virtual void diag_wet_rng(const real_t&, const real_t&)                   { assert(false); }
      virtual void diag_kappa_rng(const real_t&, const real_t&)                 { assert(false); }
      // The 3 following functions are for consecutive selection of SDs.
      // It allows the user to select SDs based on multiple characteristics, e.g. wet radius (0.5, 1) and kappa (0.1, 0.2):
      // diag_wet_rng(0.5, 1); diag_kappa_rng_cons(0.1, 0.2);
      // NOTE: the call with "cons" needs to be right after the previous call to diag_X_rng!
      //       otherwise some other function used in between can overwrite n_filtered used for selection of moments
      // NOTE2: We cannot implement this as an argument to diag_X_rng, because we would like it to default to null 
      //        and Boost Python does not work well with virtual member functions that have default arguments
      virtual void diag_dry_rng_cons(const real_t&, const real_t&)              { assert(false); }
      virtual void diag_wet_rng_cons(const real_t&, const real_t&)              { assert(false); }
      virtual void diag_kappa_rng_cons(const real_t&, const real_t&)            { assert(false); }

      virtual void diag_dry_mom(const int&)                                     { assert(false); }
      virtual void diag_wet_mom(const int&)                                     { assert(false); }
      virtual void diag_wet_mass_dens(const real_t&, const real_t&)             { assert(false); }
      virtual void diag_chem(const enum common::chem::chem_species_t&)          { assert(false); }
      virtual void diag_precip_rate()                                           { assert(false); }
      virtual void diag_kappa_mom(const int&)                                   { assert(false); }
      virtual void diag_up_mom(const int&)                                      { assert(false); }
      virtual void diag_vp_mom(const int&)                                      { assert(false); }
      virtual void diag_wp_mom(const int&)                                      { assert(false); }
      virtual std::pair<real_t, real_t> diag_up_minmax()                        { assert(false); }
      virtual std::pair<real_t, real_t> diag_vp_minmax()                        { assert(false); }
      virtual std::pair<real_t, real_t> diag_wp_minmax()                        { assert(false); }
      virtual void diag_incloud_time_mom(const int&)                            { assert(false); } // requires opts_init.diag_incloud_time==true
      virtual void diag_max_rw()                                                { assert(false); }
      virtual void diag_vel_div()                                               { assert(false); }
      virtual real_t diag_pair_separation_mean()                                { assert(false); }
      virtual int  diag_sstp_coal()                                             { assert(false); }
      virtual std::map<libcloudphxx::common::output_t, real_t> diag_puddle()    { assert(false); return std::map<libcloudphxx::common::output_t, real_t>(); }
      virtual real_t *outbuf()                                                  { assert(false); return NULL; }

      virtual void remove_wet_rng(const real_t&, const real_t&)                 { assert(false); }

      // storing a pointer to opts_init (e.g. for interrogatin about
      // dimensions in Python bindings)
      opts_init_t<real_t> *opts_init;

      // virtual destructor
      virtual ~particles_proto_t() {};

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
        const arrinfo_t<real_t> p,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_y, 
        const arrinfo_t<real_t> courant_z,
        const std::map<enum common::chem::chem_species_t, const arrinfo_t<real_t> > ambient_chem
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
        const arrinfo_t<real_t> diss_rate,
        std::map<enum common::chem::chem_species_t, arrinfo_t<real_t> > ambient_chem
      );

      void sync_in(
        arrinfo_t<real_t> th,
        arrinfo_t<real_t> rv,
        const arrinfo_t<real_t> rhod,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_y,
        const arrinfo_t<real_t> courant_z,
        const arrinfo_t<real_t> diss_rate,
        std::map<enum common::chem::chem_species_t, arrinfo_t<real_t> > ambient_chem
      );

      void step_cond(
        const opts_t<real_t> &,
        arrinfo_t<real_t> th,
        arrinfo_t<real_t> rv,
        std::map<enum common::chem::chem_species_t, arrinfo_t<real_t> > ambient_chem = std::map<enum common::chem::chem_species_t, arrinfo_t<real_t> >()
      );

      void step_async(
        const opts_t<real_t> &
      );

      // diagnostic methods
      void diag_sd_conc();
      void diag_pressure();
      void diag_temperature();
      void diag_RH();
      void diag_dry_rng(const real_t &r_mi, const real_t &r_mx);
      void diag_wet_rng(const real_t &r_mi, const real_t &r_mx);
      void diag_kappa_rng(const real_t &r_mi, const real_t &r_mx);
      void diag_dry_rng_cons(const real_t &r_mi, const real_t &r_mx);
      void diag_wet_rng_cons(const real_t &r_mi, const real_t &r_mx);
      void diag_kappa_rng_cons(const real_t &r_mi, const real_t &r_mx);
      void diag_dry_mom(const int &k);
      void diag_wet_mom(const int &k);
      void diag_kappa_mom(const int &k);
      void diag_up_mom(const int&);
      void diag_vp_mom(const int&);
      void diag_wp_mom(const int&);
      std::pair<real_t, real_t> diag_up_minmax();
      std::pair<real_t, real_t> diag_vp_minmax();
      std::pair<real_t, real_t> diag_wp_minmax();
      void diag_incloud_time_mom(const int &k);
      real_t diag_pair_separation_mean();
      void diag_wet_mass_dens(const real_t&, const real_t&);

      void diag_chem(const enum common::chem::chem_species_t&);
      void diag_rw_ge_rc();
      void diag_RH_ge_Sc();
      void diag_all();
      void diag_precip_rate();
      void diag_max_rw();
      void diag_vel_div();
      int  diag_sstp_coal();
      std::map<libcloudphxx::common::output_t, real_t> diag_puddle();
      real_t *outbuf();

      void remove_wet_rng(const real_t&, const real_t&);

      struct impl;
      std::unique_ptr<impl> pimpl;

      // constructor
      particles_t(opts_init_t<real_t> opts_init, int n_x_tot = 0); // NOTE: python bindings fail if more than one default value...

      // declare destructor to delay it's definition until impl is defined
      ~particles_t();

      // helper typedef
      typedef particles_proto_t<real_t> parent_t;
    };


    // specialization for the multi_GPU backend
    // the interface is the same as for other backends (above)
    template <typename real_t>
    struct particles_t<real_t, multi_CUDA>: particles_proto_t<real_t>
    {
      // initialisation 
      void init(
        const arrinfo_t<real_t> th,
        const arrinfo_t<real_t> rv,
        const arrinfo_t<real_t> rhod,
        const arrinfo_t<real_t> p,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_y, 
        const arrinfo_t<real_t> courant_z,
        const std::map<enum common::chem::chem_species_t, const arrinfo_t<real_t> > ambient_chem 
      );

      // time-stepping methods
      void step_sync(
        const opts_t<real_t> &,
        arrinfo_t<real_t> th,
        arrinfo_t<real_t> rv,
        const arrinfo_t<real_t> rhod     ,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_y,
        const arrinfo_t<real_t> courant_z,
        const arrinfo_t<real_t> diss_rate,
        std::map<enum common::chem::chem_species_t, arrinfo_t<real_t> > ambient_chem 
      );

      void sync_in(
        arrinfo_t<real_t> th,
        arrinfo_t<real_t> rv,
        const arrinfo_t<real_t> rhod     ,
        const arrinfo_t<real_t> courant_x,
        const arrinfo_t<real_t> courant_y,
        const arrinfo_t<real_t> courant_z,
        const arrinfo_t<real_t> diss_rate,
        std::map<enum common::chem::chem_species_t, arrinfo_t<real_t> > ambient_chem 
      );

      void step_cond(
        const opts_t<real_t> &,
        arrinfo_t<real_t> th,
        arrinfo_t<real_t> rv,
        std::map<enum common::chem::chem_species_t, arrinfo_t<real_t> > ambient_chem = std::map<enum common::chem::chem_species_t, arrinfo_t<real_t> >()
      );

      void step_async(
        const opts_t<real_t> &
      );

      // diagnostic methods
      void diag_sd_conc();
      void diag_pressure();
      void diag_temperature();
      void diag_RH();
      void diag_dry_rng(const real_t &r_mi, const real_t &r_mx);
      void diag_wet_rng(const real_t &r_mi, const real_t &r_mx);
      void diag_kappa_rng(const real_t &r_mi, const real_t &r_mx);
      void diag_dry_rng_cons(const real_t &r_mi, const real_t &r_mx);
      void diag_wet_rng_cons(const real_t &r_mi, const real_t &r_mx);
      void diag_kappa_rng_cons(const real_t &r_mi, const real_t &r_mx);
      void diag_dry_mom(const int &k);
      void diag_wet_mom(const int &k);
      void diag_kappa_mom(const int&);
      void diag_up_mom(const int&);
      void diag_vp_mom(const int&);
      void diag_wp_mom(const int&);
      std::pair<real_t, real_t> diag_up_minmax();
      std::pair<real_t, real_t> diag_vp_minmax();
      std::pair<real_t, real_t> diag_wp_minmax();
      void diag_incloud_time_mom(const int&);
      real_t diag_pair_separation_mean();
      void diag_wet_mass_dens(const real_t&, const real_t&);
      void diag_chem(const enum common::chem::chem_species_t&);
      void diag_rw_ge_rc();
      void diag_RH_ge_Sc();
      void diag_all();
      void diag_precip_rate();
      void diag_max_rw();
      void diag_vel_div();
      int  diag_sstp_coal();
      std::map<libcloudphxx::common::output_t, real_t> diag_puddle();
      real_t *outbuf();

      void remove_wet_rng(const real_t&, const real_t&);

      struct impl;
      std::unique_ptr<impl> pimpl;

      // constructors
      particles_t(opts_init_t<real_t> opts_init);

      // dtor
      ~particles_t();

      // helper typedef
      typedef particles_proto_t<real_t> parent_t;
    };
  };
};
