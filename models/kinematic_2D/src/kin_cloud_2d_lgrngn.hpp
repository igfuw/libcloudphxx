#pragma once

#include "kin_cloud_2d_common.hpp"
#include "outmom.hpp"

#include <libcloudph++/lgrngn/factory.hpp>

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

  public:
  using ix = typename ct_params_t::ix;
  using real_t = typename ct_params_t::real_t;

  protected:
  using parent_t = kin_cloud_2d_common<ct_params_t>; 

  // member fields
  std::unique_ptr<libcloudphxx::lgrngn::particles_proto_t<real_t>> prtcls;

  bool coal, sedi;

  // helper methods
  void diag()
  {
    assert(this->rank == 0);

    // recording super-droplet concentration per grid cell 
    prtcls->diag_all();
    prtcls->diag_sd_conc();
    this->record_aux("sd_conc", prtcls->outbuf());

    // recording what fell out
    auto puddle = prtcls->diag_puddle();
    for(auto elem : puddle)
    {   
       this->f_puddle << elem.first << " " << elem.second << "\n";
    }   
    this->f_puddle << "\n";

   
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
      // rw3(rd)
      int rng_num = 0;
      for (auto &rng_moms : params.out_dry)
      {
        auto &rng(rng_moms.first);
        prtcls->diag_dry_rng(rng.first / si::metres, rng.second / si::metres);
        prtcls->diag_wet_mom(3);
        this->record_aux(aux_name("rw3ofrd", rng_num, 3), prtcls->outbuf());
        rng_num++;
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
      params.cloudph_opts.coal = val ? coal : false;
      params.cloudph_opts.sedi = val ? sedi : false;
      params.cloudph_opts.RH_max = val ? 44 : 1.01; // 1% limit during spinup // TODO: specify it somewhere else, dup in blk_2m
  };

  // deals with initial supersaturation
  void hook_ante_loop(int nt)
  {
    coal = params.cloudph_opts.coal;
    sedi = params.cloudph_opts.sedi;

    parent_t::hook_ante_loop(nt); 

    // TODO: barrier?
    if (this->rank == 0) 
    {
      //save all setup parameters as attributes of setup dataset in 'const.h5' file
      this->record_aux_const("th_0",     this->setup.th_0 / si::kelvins);
      this->record_aux_const("rv_0",     this->setup.rv_0);
      this->record_aux_const("p_0",      this->setup.p_0 / si::pascals);
      this->record_aux_const("w_max",    this->setup.w_max * si::seconds / si::metres);
      this->record_aux_const("z_0",      this->setup.z_0 / si::metres);
      this->record_aux_const("X",        this->setup.X / si::metres);
      this->record_aux_const("Z",        this->setup.Z / si::metres);
      this->record_aux_const("dt",       this->setup.dt / si::seconds);
      this->record_aux_const("mean_rd1", this->setup.mean_rd1 / si::metres);
      this->record_aux_const("mean_rd2", this->setup.mean_rd2 / si::metres);
      this->record_aux_const("sdev_rd1", this->setup.sdev_rd1);
      this->record_aux_const("sdev_rd2", this->setup.sdev_rd2);
      this->record_aux_const("n1_stp",   this->setup.n1_stp * si::cubic_metres);
      this->record_aux_const("n2_stp",   this->setup.n2_stp * si::cubic_metres);
      this->record_aux_const("kappa",    this->setup.kappa);

      assert(params.backend != -1);
      assert(params.dt != 0); 

      // async does not make sense without CUDA
      if (params.backend != libcloudphxx::lgrngn::CUDA  && params.backend!= libcloudphxx::lgrngn::multi_CUDA) params.async = false;

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
        using libmpdataxx::arakawa_c::operator^;
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

	prtcls->init(
	  make_arrinfo(this->mem->advectee(ix::th)),
	  make_arrinfo(this->mem->advectee(ix::rv)),
	  make_arrinfo(this->mem->g_factor()),
          libcloudphxx::lgrngn::arrinfo_t<real_t>(),
	  make_arrinfo(Cx),
          libcloudphxx::lgrngn::arrinfo_t<real_t>(),
	  make_arrinfo(Cz)
	); 
      }

      // writing diagnostic data for the initial condition
      diag();
    }
    // TODO: barrier?
  }

#if defined(STD_FUTURE_WORKS)
  std::future<void> ftr;
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
        //!((this->timestep-1) == 0 || ((this->timestep-1) % this->outfreq == 0 && (this->timestep-1) >= this->spinup)) // .. and not after diag
      ) {
        assert(ftr.valid());
        ftr.get();
      } else assert(!ftr.valid());
#endif

      // running synchronous stuff
      prtcls->step_sync(
        params.cloudph_opts,
        make_arrinfo(this->mem->advectee(ix::th)),
        make_arrinfo(this->mem->advectee(ix::rv))
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
              dynamic_cast<particles_t<real_t, multi_CUDA>*>(prtcls.get()),
              params.cloudph_opts
            );
          else if(params.backend == CUDA)
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
      //if (this->timestep == 0 || (this->timestep % this->outfreq == 0 && this->timestep >= this->spinup))
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
    outmom_t<real_t> out_dry, out_wet;
  };

  protected:

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
