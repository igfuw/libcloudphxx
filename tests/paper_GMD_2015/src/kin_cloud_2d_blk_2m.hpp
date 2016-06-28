#include "kin_cloud_2d_common.hpp"

#include <libcloudph++/blk_2m/options.hpp>
#include <libcloudph++/blk_2m/rhs_cellwise.hpp>
#include <libcloudph++/blk_2m/rhs_columnwise.hpp>

template <class ct_params_t>
class kin_cloud_2d_blk_2m : public kin_cloud_2d_common<ct_params_t>
{
  using parent_t = kin_cloud_2d_common<ct_params_t>;

  public:
  using ix = typename ct_params_t::ix;
  using real_t = typename ct_params_t::real_t;
  private:

  void zero_if_uninitialised(int e)  //TODO move to common
  {
    if (!std::isfinite(sum(this->state(e)(this->ijk)))) 
    this->state(e)(this->ijk) = 0;
  }

  void update_rhs(
    libmpdataxx::arrvec_t<typename parent_t::arr_t> &rhs,
    const typename parent_t::real_t &dt,
    const int &at 
  ) {
    parent_t::update_rhs(rhs, dt, at);

    this->mem->barrier(); // TODO: if neccesarry, then move to adv_rhs/....hpp

    // cell-wise
    {
      auto
	dot_th = rhs.at(ix::th)(this->ijk),
	dot_rv = rhs.at(ix::rv)(this->ijk),
	dot_rc = rhs.at(ix::rc)(this->ijk),
	dot_rr = rhs.at(ix::rr)(this->ijk),
	dot_nc = rhs.at(ix::nc)(this->ijk),
	dot_nr = rhs.at(ix::nr)(this->ijk);
      const auto
        rhod   = (*this->mem->G)(this->ijk),
        th     = this->state(ix::th)(this->ijk),
        rv     = this->state(ix::rv)(this->ijk),
        rc     = this->state(ix::rc)(this->ijk),
        rr     = this->state(ix::rr)(this->ijk),
        nc     = this->state(ix::nc)(this->ijk),
        nr     = this->state(ix::nr)(this->ijk);

      libcloudphxx::blk_2m::rhs_cellwise<real_t>(
        opts, dot_th, dot_rv, dot_rc, dot_nc, dot_rr, dot_nr,
	rhod,     th,     rv,     rc,     nc,     rr,     nr,
        this->dt
      );
    }

    // column-wise
    for (int i = this->i.first(); i <= this->i.last(); ++i)
    {
      auto
	dot_rr = rhs.at(ix::rr)(i, this->j),
	dot_nr = rhs.at(ix::nr)(i, this->j);
      const auto
        rhod   = (*this->mem->G)(i, this->j),
	rr     = this->state(ix::rr)(i, this->j),
	nr     = this->state(ix::nr)(i, this->j);

      libcloudphxx::blk_2m::rhs_columnwise<real_t>(
        opts, dot_rr, dot_nr, 
        rhod,     rr,     nr,  
	this->dt,
	this->dz
      );
    }

    this->mem->barrier(); // TODO: if needed, move to adv+rhs
  }

  libcloudphxx::blk_2m::opts_t<real_t> opts;

  protected:

  void hook_ante_loop(int nt)
  {
    // if uninitialised fill with zeros
    zero_if_uninitialised(ix::rc);
    zero_if_uninitialised(ix::nc);
    zero_if_uninitialised(ix::rr);
    zero_if_uninitialised(ix::nr);

    parent_t::hook_ante_loop(nt); 
  }

  // spinup stuff
  bool get_rain() { return opts.acnv; }
  void set_rain(bool val) 
  { 
    opts.acnv = val; 
    opts.RH_max = val ? 44 : 1.01; // 1% limit during spinup
  };

  public:

  struct rt_params_t : parent_t::rt_params_t 
  { 
    libcloudphxx::blk_2m::opts_t<real_t> cloudph_opts;
  };

  // ctor
  kin_cloud_2d_blk_2m( 
    typename parent_t::ctor_args_t args, 
    const rt_params_t &p
  ) : 
    parent_t(args, p),
    opts(p.cloudph_opts)
  { 
    assert(p.dt != 0);
  }  
};

