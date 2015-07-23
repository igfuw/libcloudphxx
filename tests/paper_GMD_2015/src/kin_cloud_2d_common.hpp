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
  bool relax_th_rv;

  // spinup stuff
  virtual bool get_rain() = 0;
  virtual void set_rain(bool) = 0;

  void hook_ante_loop(int nt) 
  {
    if (get_rain() == false) spinup = 0; // spinup does not make sense without autoconversion  (TODO: issue a warning?)
    if (spinup > 0) set_rain(false);

    parent_t::hook_ante_loop(nt); 
  }

  void hook_ante_step()
  {
    // turn autoconversion on only after spinup (if spinup was specified)
    if (spinup != 0 && spinup == this->timestep) set_rain(true);

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

    // relaxation terms
    if(relax_th_rv)
    {
      // computed level-wise
      for (int j = this->j.first(); j <= this->j.last(); ++j)
      {  
        const auto tau = icmw8_case1::tau_rlx / si::seconds * exp(j * dz / icmw8_case1::z_rlx * si::metres);

        for(auto a: std::list<int>({ix::th, ix::rv}))
        {
          const auto &psi = this->state(a);
          const auto psi_mean = this->mem->sum(psi, this->i, rng_t(j, j), false)  /  this->mem->grid_size[0];
          rhs.at(a)(this->i, j) = - (psi(this->i, j) - psi_mean) / tau;
        }
      }
    }
  }

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    typename ct_params_t::real_t dx = 0, dz = 0;
    int spinup = 0; // number of timesteps during which autoconversion is to be turned off
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
    relax_th_rv(p.relax_th_rv)
  {
    assert(dx != 0);
    assert(dz != 0);
  }  
};
