// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::adjust_timesteps(const real_t &_dt)
    {
      if(_dt > 0 && !opts_init.variable_dt_switch) throw std::runtime_error("opts.dt specified, but opts_init.variable_dt_switch is false.");
      // dt defined in opts_init can be overriden by dt passed to this function
      dt = _dt > 0 ? _dt : opts_init.dt;
      // then, number of substeps is adjusted to get desired dt for specific processes
      sstp_cond = _dt > 0 ? std::ceil(opts_init.sstp_cond * _dt / opts_init.dt) : opts_init.sstp_cond;
      sstp_coal = _dt > 0 ? std::ceil(opts_init.sstp_coal * _dt / opts_init.dt) : opts_init.sstp_coal;
      sstp_chem = _dt > 0 ? std::ceil(opts_init.sstp_chem * _dt / opts_init.dt) : opts_init.sstp_chem;
    }
  };
};
