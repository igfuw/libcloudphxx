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
      if(_dt > 0 && !opts_init.variable_dt_switch) throw std::runtime_error("libcloudph++: opts.dt specified, but opts_init.variable_dt_switch is false.");
      // then, number of substeps is adjusted to get desired dt for specific processes, but only if initially they are > 1
      sstp_cond = _dt > 0     && sstp_cond > 1     ? std::ceil(sstp_cond * _dt / dt)     : sstp_cond;
      sstp_coal = _dt > 0     && sstp_coal > 1     ? std::ceil(sstp_coal * _dt / dt)     : sstp_coal;
      sstp_chem = _dt > 0     && sstp_chem > 1     ? std::ceil(sstp_chem * _dt / dt)     : sstp_chem;
      sstp_cond_act = _dt > 0 && sstp_cond_act > 1 ? std::ceil(sstp_cond_act * _dt / dt) : sstp_cond_act;
      // dt defined in opts_init can be overriden by dt passed to this function
      dt = _dt > 0 ? _dt : opts_init.dt;
    }
  };
};
