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
    // calc the SGS turbulent velocity component
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_sgs_vel(const real_t &dt, const bool only_vertical)
    {   
      switch(opts_init.sgs_adve)
      {
        case sgs_adve_t::GA17:
          hskpng_sgs_vel_GA17(dt, only_vertical);
          break;
        case sgs_adve_t::ST_periodic:
          hskpng_sgs_vel_ST_periodic(dt, only_vertical);
          break;
        default:
          throw std::runtime_error("Unreckognized value of opts_init.sgs_adve.");
      }
    }
  };
};
