#pragma once

#include "icmw8_case1.hpp"
#include "kin_cloud_2d_lgrngn.hpp"
#include "outmom.hpp"

#include <libcloudph++/lgrngn/factory.hpp>
#include <libcloudph++/lgrngn/chem.hpp>

// @brief a minimalistic kinematic cloud model with lagrangian microphysics
//        built on top of the mpdata_2d solver (by extending it with
//        custom hook_ante_loop() and hook_post_step() methods)
template <class ct_params_t>
class kin_cloud_2d_lgrngn_chem : public kin_cloud_2d_lgrngn<ct_params_t>
{
  // note: lgrngn has no rhs terms - just adjustments (but there might be extrinsic rhs terms)
  using parent_t = kin_cloud_2d_lgrngn<ct_params_t>; 

  public:
  using ix = typename ct_params_t::ix;
  using real_t = typename ct_params_t::real_t;
  private:

  // recording mass of H and S_VI in wet radius bins
  void diag_pH()
  {
    assert(this->rank == 0);
    {
      int rng_num = 0;
      for (auto &rng_moms : params.out_wet_pH)
      {
        //wet
        auto &rng(rng_moms.first);
        parent_t::prtcls->diag_wet_rng(rng.first / si::metres, rng.second / si::metres);
        for (auto &mom : rng_moms.second)
        {
          parent_t::prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::S_VI);
          this->record_aux(this->aux_name("chem_S_VI_rw", rng_num, mom), parent_t::prtcls->outbuf());
 
          parent_t::prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::H);
          this->record_aux(this->aux_name("chem_H_rw", rng_num, mom), parent_t::prtcls->outbuf());
        }
        rng_num++;
      }
    }
  }

  void diag_chem()
  {
    assert(this->rank == 0);
    {
      // chem
      for (auto &rng_moms : params.out_chem)
      {
        auto &rng(rng_moms.first);
        parent_t::prtcls->diag_dry_rng(rng.first / si::metres, rng.second / si::metres);
         
        //TODO: for (auto &chem_enum : libcloudphxx::lgrngn::chem_species_t)
        {
          parent_t::prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::SO2);
          this->record_aux("chem_S_IV_aq", parent_t::prtcls->outbuf());
          parent_t::prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::S_VI);
          this->record_aux("chem_S_VI_aq", parent_t::prtcls->outbuf());

          parent_t::prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::O3);
          this->record_aux("chem_O3_aq", parent_t::prtcls->outbuf());
          parent_t::prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::H2O2);
          this->record_aux("chem_H2O2_aq", parent_t::prtcls->outbuf());

          parent_t::prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::H);
          this->record_aux("chem_H_aq", parent_t::prtcls->outbuf());

          parent_t::prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::CO2);
          this->record_aux("chem_C_IV_aq", parent_t::prtcls->outbuf());

          parent_t::prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::NH3);
          this->record_aux("chem_N_III_aq", parent_t::prtcls->outbuf());

          parent_t::prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::HNO3);
          this->record_aux("chem_N_V_aq", parent_t::prtcls->outbuf());
        }
      }
    }
  }
  protected:

  void set_chem(bool val) 
  {
    if (val && params.cloudph_opts_init.chem_switch == true) 
      params.cloudph_opts.chem_rct = val;
  };

  void set_rain(bool val)
  {
    parent_t::set_rain(val);
    set_chem(val);
  }

  // initial supersaturation + init for chem 
  void hook_ante_loop(int nt)
  {
    bool coal, sedi;
    coal = params.cloudph_opts.coal;
    sedi = params.cloudph_opts.sedi;

    {
      blitz::secondIndex j;
      namespace molar_mass  = libcloudphxx::common::molar_mass;

      // initialise trace gases profiles
      this->state(ix::SO2g)  = 
        config::mixr_helper(this->setup)(j * params.dz) 
        * (this->setup.SO2_g_0  * molar_mass::M_SO2<real_t>()  * si::moles / si::kilograms);
      this->state(ix::O3g)   = 
        config::mixr_helper(this->setup)(j * params.dz) 
        * (this->setup.O3_g_0   * molar_mass::M_O3<real_t>()   * si::moles / si::kilograms);
      this->state(ix::H2O2g) = 
        config::mixr_helper(this->setup)(j * params.dz) 
        * (this->setup.H2O2_g_0 * molar_mass::M_H2O2<real_t>() * si::moles / si::kilograms);
      this->state(ix::CO2g)  = 
        config::mixr_helper(this->setup)(j * params.dz) 
        * (this->setup.CO2_g_0  * molar_mass::M_CO2<real_t>()  * si::moles / si::kilograms);
      this->state(ix::NH3g)  = 
        config::mixr_helper(this->setup)(j * params.dz) 
        * (this->setup.NH3_g_0  * molar_mass::M_NH3<real_t>()  * si::moles / si::kilograms);
      this->state(ix::HNO3g) = 
        config::mixr_helper(this->setup)(j * params.dz) 
        * (this->setup.HNO3_g_0 * molar_mass::M_HNO3<real_t>() * si::moles / si::kilograms);
    }

    parent_t::parent_t::hook_ante_loop(nt); 

   // TODO: barrier?
   // TODO - make a lagrangian_common class for particles with and without chem 
    if (this->rank == 0) 
    {
      assert(params.backend != -1);
      assert(params.dt != 0); 

      // async does not make sense without CUDA
      if (params.backend != libcloudphxx::lgrngn::CUDA) params.async = false;

      params.cloudph_opts_init.dt = params.dt; // advection timestep = microphysics timestep
      params.cloudph_opts_init.dx = params.dx;
      params.cloudph_opts_init.dz = params.dz;


      // libmpdata++'s grid interpretation
      params.cloudph_opts_init.x0 = params.dx / 2;
      params.cloudph_opts_init.z0 = params.dz / 2;
      params.cloudph_opts_init.x1 = (this->mem->grid_size[0].length() - .5) * params.dx;
      params.cloudph_opts_init.z1 = (this->mem->grid_size[1].length() - .5) * params.dz;

      parent_t::prtcls.reset(libcloudphxx::lgrngn::factory<real_t>(
        (libcloudphxx::lgrngn::backend_t)params.backend, 
        params.cloudph_opts_init
      ));

      {
        using libmpdataxx::arakawa_c::h;
        // temporarily Cx & Cz are multiplied by rhod ...
        auto 
          Cx = this->mem->GC[0](
            this->mem->grid_size[0]^h, 
            this->mem->grid_size[1]
          ).reindex({0,0}).copy(),
          Cz = this->mem->GC[1](
            this->mem->grid_size[0], 
            this->mem->grid_size[1]^h
          ).reindex({0,0}).copy();

        // ... and now dividing them by rhod (z=0 is located at j=1/2)
        {
          blitz::secondIndex j;
          Cx /= config::rhod(this->setup)(   j     * this->dj);
          Cz /= config::rhod(this->setup)((j - .5) * this->dj);
        }

        assert(params.cloudph_opts_init.chem_switch == true);

        std::map<enum libcloudphxx::lgrngn::chem_species_t, const libcloudphxx::lgrngn::arrinfo_t<real_t> > ambient_chem_init;
        boost::assign::insert(ambient_chem_init)
          (libcloudphxx::lgrngn::chem_species_t::SO2,  this->make_arrinfo(this->mem->advectee(ix::SO2g)))
          (libcloudphxx::lgrngn::chem_species_t::O3,   this->make_arrinfo(this->mem->advectee(ix::O3g)))
          (libcloudphxx::lgrngn::chem_species_t::H2O2, this->make_arrinfo(this->mem->advectee(ix::H2O2g)))
          (libcloudphxx::lgrngn::chem_species_t::CO2,  this->make_arrinfo(this->mem->advectee(ix::CO2g)))
          (libcloudphxx::lgrngn::chem_species_t::NH3,  this->make_arrinfo(this->mem->advectee(ix::NH3g)))
          (libcloudphxx::lgrngn::chem_species_t::HNO3, this->make_arrinfo(this->mem->advectee(ix::HNO3g)));

        // TODO!!! 
        // why we need those is beyond me, but apparently without it
        // strides in ambient_chem are wrong for some chemical species
        // the problem is not present in Debug mode
        libcloudphxx::lgrngn::arrinfo_t<real_t> th(this->make_arrinfo(this->mem->advectee(ix::th)));
        libcloudphxx::lgrngn::arrinfo_t<real_t> rv(this->make_arrinfo(this->mem->advectee(ix::rv)));
        libcloudphxx::lgrngn::arrinfo_t<real_t> g( this->make_arrinfo(this->mem->g_factor()));

	parent_t::prtcls->init(
	  this->make_arrinfo(this->mem->advectee(ix::th)),
	  this->make_arrinfo(this->mem->advectee(ix::rv)),
	  this->make_arrinfo(this->mem->g_factor()),
	  this->make_arrinfo(Cx),
          libcloudphxx::lgrngn::arrinfo_t<real_t>(),
	  this->make_arrinfo(Cz),
          ambient_chem_init
	); 
      }

      // writing diagnostic data for the initial condition

      parent_t::params.out_wet = params.out_wet;
      parent_t::params.out_dry = params.out_dry;

      parent_t::diag();
      diag_chem();
      diag_pH();
    }
    // TODO: barrier?
  }

#if defined(STD_FUTURE_WORKS)
  std::future<real_t> ftr;
#endif

  // 
  void hook_post_step()
  {
    parent_t::parent_t::hook_post_step(); // includes output

    this->mem->barrier();

    if (this->rank == 0) 
    {
      // assuring previous async step finished ...
#if defined(STD_FUTURE_WORKS)
      if (
        params.async && 
        this->timestep != 0 && // ... but not in first timestep ...
        ((this->timestep - 1) % this->outfreq != 0) // ... and not after diag call
      ) {
        assert(ftr.valid());
        ftr.get();
      } else assert(!ftr.valid());
#endif

      assert(params.cloudph_opts_init.chem_switch == true);

      std::map<enum libcloudphxx::lgrngn::chem_species_t, libcloudphxx::lgrngn::arrinfo_t<real_t> > ambient_chem;
      boost::assign::insert(ambient_chem)
        (libcloudphxx::lgrngn::chem_species_t::SO2,  this->make_arrinfo(this->mem->advectee(ix::SO2g)))
        (libcloudphxx::lgrngn::chem_species_t::O3,   this->make_arrinfo(this->mem->advectee(ix::O3g)))
        (libcloudphxx::lgrngn::chem_species_t::H2O2, this->make_arrinfo(this->mem->advectee(ix::H2O2g)))
        (libcloudphxx::lgrngn::chem_species_t::CO2,  this->make_arrinfo(this->mem->advectee(ix::CO2g)))
        (libcloudphxx::lgrngn::chem_species_t::NH3,  this->make_arrinfo(this->mem->advectee(ix::NH3g)))
        (libcloudphxx::lgrngn::chem_species_t::HNO3, this->make_arrinfo(this->mem->advectee(ix::HNO3g)));

      // running synchronous stuff
      parent_t::prtcls->step_sync(
        params.cloudph_opts,
        this->make_arrinfo(this->mem->advectee(ix::th)),
        this->make_arrinfo(this->mem->advectee(ix::rv)),
        libcloudphxx::lgrngn::arrinfo_t<real_t>(),
        libcloudphxx::lgrngn::arrinfo_t<real_t>(),
        libcloudphxx::lgrngn::arrinfo_t<real_t>(),
        libcloudphxx::lgrngn::arrinfo_t<real_t>(),
        ambient_chem
      ); 

      // running asynchronous stuff
      {
        using libcloudphxx::lgrngn::particles_t;
        using libcloudphxx::lgrngn::CUDA;
        using libcloudphxx::lgrngn::multi_CUDA;

#if defined(STD_FUTURE_WORKS)
        if (params.async)
        {
          assert(!ftr.valid());

          if(params.backend == multi_CUDA)
            ftr = std::async(
              std::launch::async, 
              &particles_t<real_t, multi_CUDA>::step_async, 
              dynamic_cast<particles_t<real_t, multi_CUDA>*>(parent_t::prtcls.get()),
              params.cloudph_opts
            );
          else if(params.backend == CUDA)
            ftr = std::async(
              std::launch::async, 
              &particles_t<real_t, CUDA>::step_async, 
              dynamic_cast<particles_t<real_t, CUDA>*>(parent_t::prtcls.get()),
              params.cloudph_opts
            );
          assert(ftr.valid());
        } else 
#endif
          parent_t::prtcls->step_async(params.cloudph_opts);
      }

      // performing diagnostics
      if (this->timestep % this->outfreq == 0) 
      { 
#if defined(STD_FUTURE_WORKS)
        if (params.async)
        {
          assert(ftr.valid());
          ftr.get();
        }
#endif
        parent_t::diag();
        diag_chem();
        diag_pH();
      }
    }

    this->mem->barrier();
  }

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    int backend = -1;
    bool async = true;
    libcloudphxx::lgrngn::opts_t<real_t> cloudph_opts;
    libcloudphxx::lgrngn::opts_init_t<real_t> cloudph_opts_init;
    outmom_t<real_t> out_dry, out_wet, out_chem, out_wet_pH;
  };

  private:

  // per-thread copy of params
  rt_params_t params;

  public:

  // ctor
  kin_cloud_2d_lgrngn_chem( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) : 
    parent_t(args, p),
    params(p)
  {
    // delaying any initialisation to ante_loop as rank() does not function within ctor! // TODO: not anymore!!!
    // TODO: equip rank() in libmpdata with an assert() checking if not in serial block
  }  
};
