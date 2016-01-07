#pragma once

#include <vector>

#include "kin_cloud_2d_common.hpp"
#include "outmom.hpp"

#include <libcloudph++/lgrngn/factory.hpp>
#include <libcloudph++/lgrngn/chem.hpp>

#if defined(STD_FUTURE_WORKS)
#  include <future>
#endif

// @brief a minimalistic kinematic cloud model with lagrangian microphysics
//        built on top of the mpdata_2d solver (by extending it with
//        custom hook_ante_loop() and hook_post_step() methods)
template <class ct_params_t>
class kin_cloud_2d_lgrngn : public kin_cloud_2d_common<ct_params_t>
{
  // note: lgrngn has no rhs terms - just adjustments (but there might be extrinsic rhs terms)
  using parent_t = kin_cloud_2d_common<ct_params_t>; 

  public:
  using ix = typename ct_params_t::ix;
  using real_t = typename ct_params_t::real_t;
  private:

  // member fields
  std::unique_ptr<libcloudphxx::lgrngn::particles_proto_t<real_t>> prtcls;

  // helper methods
  void diag()
  {
    assert(this->rank == 0);

    // recording super-droplet concentration per grid cell 
    prtcls->diag_sd_conc();
    this->record_aux("sd_conc", prtcls->outbuf());
   
    // recording requested statistical moments
    {
      // dry
      int rng_num = 0;
      for (auto &rng_moms : params.out_dry)
      {
        auto &rng(rng_moms.first);
        prtcls->diag_dry_rng(rng.first / si::metres, rng.second / si::metres);
        for (auto &mom : rng_moms.second)
        {
          prtcls->diag_dry_mom(mom);
          this->record_aux(aux_name("rd", rng_num, mom), prtcls->outbuf());
        }
        rng_num++;
      }
    }
    {
      // wet
      int rng_num = 0;
      for (auto &rng_moms : params.out_wet)
      {
        auto &rng(rng_moms.first);
        prtcls->diag_wet_rng(rng.first / si::metres, rng.second / si::metres);
        for (auto &mom : rng_moms.second)
        {
          prtcls->diag_wet_mom(mom);
          this->record_aux(aux_name("rw", rng_num, mom), prtcls->outbuf());
        }
        rng_num++;
      }
    }

    {
      // chem
      for (auto &rng_moms : params.out_chem)
      {
        auto &rng(rng_moms.first);
        prtcls->diag_dry_rng(rng.first / si::metres, rng.second / si::metres);
         
        //TODO
        //for (auto &chem_enum : libcloudphxx::lgrngn::chem_species_t)
        {
          prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::SO2);
          this->record_aux("chem_SO2_aq", prtcls->outbuf());
          prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::HSO3);
          this->record_aux("chem_HSO3_aq", prtcls->outbuf());
          prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::SO3);
          this->record_aux("chem_SO3_aq", prtcls->outbuf());
          prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::HSO4);
          this->record_aux("chem_HSO4_aq", prtcls->outbuf());
          prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::SO4);
          this->record_aux("chem_SO4_aq", prtcls->outbuf());
          prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::S_VI);
          this->record_aux("chem_S_VI_aq", prtcls->outbuf());

          prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::O3);
          this->record_aux("chem_O3_aq", prtcls->outbuf());
          prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::H2O2);
          this->record_aux("chem_H2O2_aq", prtcls->outbuf());

          prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::H);
          this->record_aux("chem_H_aq", prtcls->outbuf());
          prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::OH);
          this->record_aux("chem_OH_aq", prtcls->outbuf());

          prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::CO2);
          this->record_aux("chem_CO2_aq", prtcls->outbuf());
          prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::HCO3);
          this->record_aux("chem_HCO3_aq", prtcls->outbuf());
          prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::CO3);
          this->record_aux("chem_CO3_aq", prtcls->outbuf());

          prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::NH3);
          this->record_aux("chem_NH3_aq", prtcls->outbuf());
          prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::NH4);
          this->record_aux("chem_NH4_aq", prtcls->outbuf());

          prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::HNO3);
          this->record_aux("chem_HNO3_aq", prtcls->outbuf());
          prtcls->diag_chem(libcloudphxx::lgrngn::chem_species_t::NO3);
          this->record_aux("chem_NO3_aq", prtcls->outbuf());
        }
      }
    }
  } 

  libcloudphxx::lgrngn::arrinfo_t<real_t> make_arrinfo(
    typename parent_t::arr_t arr
  ) {
    return libcloudphxx::lgrngn::arrinfo_t<real_t>(
      arr.dataZero(), 
      arr.stride().data()
    );
  }

  std::string aux_name(
    const std::string pfx, 
    const int rng,
    const int mom
  )
  { 
    std::ostringstream tmp;
    tmp << pfx << "_rng" << std::setw(3) << std::setfill('0') << rng << "_mom" << mom;
    return tmp.str();
  }

  protected:

  bool get_rain() { return params.cloudph_opts.coal && params.cloudph_opts.sedi; }
  void set_rain(bool val) 
  { 
    params.cloudph_opts.coal = params.cloudph_opts.sedi = val;
    params.cloudph_opts.RH_max = val ? 44 : 1.01; // 1% limit during spinup // TODO: specify it somewhere else, dup in blk_2m
  };
  void set_chem(bool val) 
  {
    if (val && params.cloudph_opts_init.chem_switch == true) 
      params.cloudph_opts.chem_rct = val;
  };

  // deals with initial supersaturation
  void hook_ante_loop(int nt)
  {
    parent_t::hook_ante_loop(nt); 

    // TODO: barrier?
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

      prtcls.reset(libcloudphxx::lgrngn::factory<real_t>(
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
          Cx /= chem_case::rhod()(   j     * this->dj);
          Cz /= chem_case::rhod()((j - .5) * this->dj);
        }

        std::map<enum libcloudphxx::lgrngn::chem_species_t, const libcloudphxx::lgrngn::arrinfo_t<real_t> > ambient_chem_init;
        boost::assign::insert(ambient_chem_init)
          (libcloudphxx::lgrngn::chem_species_t::SO2,  make_arrinfo(this->mem->advectee(ix::SO2g)))
          (libcloudphxx::lgrngn::chem_species_t::O3,   make_arrinfo(this->mem->advectee(ix::O3g)))
          (libcloudphxx::lgrngn::chem_species_t::H2O2, make_arrinfo(this->mem->advectee(ix::H2O2g)))
          (libcloudphxx::lgrngn::chem_species_t::CO2,  make_arrinfo(this->mem->advectee(ix::CO2g)))
          (libcloudphxx::lgrngn::chem_species_t::NH3,  make_arrinfo(this->mem->advectee(ix::NH3g)))
          (libcloudphxx::lgrngn::chem_species_t::HNO3, make_arrinfo(this->mem->advectee(ix::HNO3g)));
/*
        std::cerr << "amb hno3 strides map " << ambient_chem_init.at((libcloudphxx::lgrngn::chem_species_t)0).strides << std::endl;
        std::cerr << "amb hno3 strides[0] map " << ambient_chem_init.at((libcloudphxx::lgrngn::chem_species_t)0).strides[0] << std::endl;

        std::cerr << "amb hno3 strdes[0] "<< this->mem->advectee(ix::HNO3g).stride().data()[0] << std::endl;
        std::cerr <<  "amb hno3 strdes[1] "<<this->mem->advectee(ix::HNO3g).stride().data()[1] << std::endl;

        std::cerr << "amb so2 strdes[0] "<< this->mem->advectee(ix::SO2g).stride().data()[0] << std::endl;
        std::cerr <<  "amb so2 strdes[1] "<<this->mem->advectee(ix::SO2g).stride().data()[1] << std::endl;
*/
        // TODO!!! 
        // why we need those is beyond me, but apparently without it
        // strides in ambient_chem are wrong for some chemical species
        // the problem is not present in Debug mode
        libcloudphxx::lgrngn::arrinfo_t<real_t> th(make_arrinfo(this->mem->advectee(ix::th)));
        libcloudphxx::lgrngn::arrinfo_t<real_t> rv(make_arrinfo(this->mem->advectee(ix::rv)));
        libcloudphxx::lgrngn::arrinfo_t<real_t> g(make_arrinfo(this->mem->g_factor()));

	prtcls->init(
	  make_arrinfo(this->mem->advectee(ix::th)),
	  make_arrinfo(this->mem->advectee(ix::rv)),
	  make_arrinfo(this->mem->g_factor()),
	  make_arrinfo(Cx),
          libcloudphxx::lgrngn::arrinfo_t<real_t>(),
	  make_arrinfo(Cz),
          ambient_chem_init
	); 
      }

      // writing diagnostic data for the initial condition
      diag();
    }
    // TODO: barrier?
  }

#if defined(STD_FUTURE_WORKS)
  std::future<real_t> ftr;
#endif

  // 
  void hook_post_step()
  {
    parent_t::hook_post_step(); // includes output

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

     std::map<enum libcloudphxx::lgrngn::chem_species_t, libcloudphxx::lgrngn::arrinfo_t<real_t> > ambient_chem;
        boost::assign::insert(ambient_chem)
          (libcloudphxx::lgrngn::chem_species_t::SO2,  make_arrinfo(this->mem->advectee(ix::SO2g)))
          (libcloudphxx::lgrngn::chem_species_t::O3,   make_arrinfo(this->mem->advectee(ix::O3g)))
          (libcloudphxx::lgrngn::chem_species_t::H2O2, make_arrinfo(this->mem->advectee(ix::H2O2g)))
          (libcloudphxx::lgrngn::chem_species_t::CO2,  make_arrinfo(this->mem->advectee(ix::CO2g)))
          (libcloudphxx::lgrngn::chem_species_t::NH3,  make_arrinfo(this->mem->advectee(ix::NH3g)))
          (libcloudphxx::lgrngn::chem_species_t::HNO3, make_arrinfo(this->mem->advectee(ix::HNO3g)));

      // running synchronous stuff
      prtcls->step_sync(
        params.cloudph_opts,
        make_arrinfo(this->mem->advectee(ix::th)),
        make_arrinfo(this->mem->advectee(ix::rv)),
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

#if defined(STD_FUTURE_WORKS)
        if (params.async)
        {
          assert(!ftr.valid());
          ftr = std::async(
            std::launch::async, 
            &particles_t<real_t, CUDA>::step_async, 
            dynamic_cast<particles_t<real_t, CUDA>*>(prtcls.get()),
            params.cloudph_opts
          );
          assert(ftr.valid());
        } else 
#endif
          prtcls->step_async(params.cloudph_opts);
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
        diag();
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
    outmom_t<real_t> out_dry, out_wet, out_chem;
  };

  private:

  // per-thread copy of params
  rt_params_t params;

  public:

  // ctor
  kin_cloud_2d_lgrngn( 
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
