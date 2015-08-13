#pragma once

#include <libmpdata++/solvers/mpdata_rhs.hpp>
#include <libmpdata++/output/hdf5.hpp>

using namespace libmpdataxx; // TODO: get rid of it?

template <class ct_params_t>
class kin_cloud_2d_common : public 
  output::hdf5<
    solvers::mpdata_rhs<ct_params_t>
  >
{
  using parent_t = output::hdf5<solvers::mpdata_rhs<ct_params_t>>;

  protected:

  typename ct_params_t::real_t dx, dz; // 0->dx, 1->dy ! TODO
  int spinup; // number of timesteps
  int src_t0, src_t1; //time interval of source

  // relaxation stuff
  bool relax_th_rv;
  blitz::Array<typename ct_params_t::real_t, 1> th_eq, rv_eq;
  
  // spinup stuff
  virtual bool get_rain() = 0;
  virtual void set_rain(bool) = 0;

  // source stuff
  virtual bool get_src() = 0;
  virtual void set_src(bool) = 0;

  void hook_ante_loop(int nt) 
  {
    if (get_rain() == false) 
    {
      // spinup and relaxation do not make sense without autoconversion  (TODO: issue a warning?)
      spinup = relax_th_rv = 0;      
    }
    if (spinup > 0) set_rain(false);
    if (src_t0 > 0) set_src(false);

    parent_t::hook_ante_loop(nt); 
  }

  void hook_ante_step()
  {
    if (spinup != 0 && spinup == this->timestep)
    {
      // turn autoconversion on only after spinup (if spinup was specified)
      set_rain(true);
    }

    if (src_t0 != 0 && src_t0 == this->timestep)
    {
      set_src(true);
    }
    if (src_t1 != 0 && src_t1 == this->timestep)
    {
      set_src(false);
    }

    using ix = typename ct_params_t::ix;
    if(relax_th_rv && spinup == this->timestep)
    {
      // save horizontal means of th and rv after spinup
      // they will be the relaxation goals
      for (int j = this->j.first(); j <= this->j.last(); ++j)
      {  
        th_eq(j) = this->mem->sum(this->state(ix::th), rng_t(this->i.first()+1, this->i.last()), rng_t(j, j), false)  /  (this->mem->grid_size[0].length() - 1);
        rv_eq(j) = this->mem->sum(this->state(ix::rv), rng_t(this->i.first()+1, this->i.last()), rng_t(j, j), false)  /  (this->mem->grid_size[0].length() - 1);
      }
    }

    parent_t::hook_ante_step(); 
  }


  void update_rhs(
    arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at 
  )   
  {   
    parent_t::update_rhs(rhs, dt, at);
    using ix = typename ct_params_t::ix;

    // relaxation terms; added only after spinup, when get_rain returns true
    if(relax_th_rv && get_rain())
    {
      // computed level-wise
      for (int j = this->j.first(); j <= this->j.last(); ++j)
      {  
        const auto tau = icmw8_case1::tau_rlx / si::seconds * exp(j * dz / icmw8_case1::z_rlx * si::metres);

        for(auto a: std::list<int>({ix::th, ix::rv}))
        {
          const auto &psi = this->state(a);
          // relax horizontal mean
          /*
          const auto psi_mean = this->mem->sum(psi, rng_t(this->i.first()+1, this->i.last()), rng_t(j, j), false)  /  (this->mem->grid_size[0].length() - 1);
          if(a == ix::th)
            rhs.at(a)(this->i, j) =  (th_eq(j) - psi_mean) / tau;
          else
            rhs.at(a)(this->i, j) =  (rv_eq(j) - psi_mean) / tau;
          */
          // relax each cell 
          if(a == ix::th)
            rhs.at(a)(this->i, j) =  (th_eq(j) - psi(this->i, j)) / tau;
          else
            rhs.at(a)(this->i, j) =  (rv_eq(j) - psi(this->i, j)) / tau;
        }
      }
    }
  }

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    typename ct_params_t::real_t dx = 0, dz = 0;
    int spinup = 0; // number of timesteps during which autoconversion is to be turned off
    int src_t0, src_t1;
    bool relax_th_rv = true;
  };

  // ctor
  kin_cloud_2d_common( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) : 
    parent_t(args, p),
    dx(p.dx),
    dz(p.dz),
    spinup(p.spinup),
    src_t0(p.src_t0),
    src_t1(p.src_t1),
    relax_th_rv(p.relax_th_rv),
    th_eq(this->mem->grid_size[1].length()),
    rv_eq(this->mem->grid_size[1].length())
  {
    assert(dx != 0);
    assert(dz != 0);
  }  
};
