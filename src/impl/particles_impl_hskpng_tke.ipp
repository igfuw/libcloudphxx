// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <libcloudph++/common/GA17_turbulence.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template<class real_t>
      struct common__turbulence__tke
      {
        const quantity<si::length, real_t> lambda;
        common__turbulence__tke(const real_t &lambda):
          lambda(lambda * si::metres){}

        BOOST_GPU_ENABLED
        real_t operator()(const real_t &diss_rate)
        {
          return common::GA17_turbulence::tke(
            diss_rate * si::metres * si::metres / si::seconds / si::seconds / si::seconds,
            lambda) / si::metres / si::metres * si::seconds * si::seconds;
        }
      };
    };
    // calc the SGS TKE, in place of dissipation rate
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_tke()
    {   
      thrust::transform(diss_rate.begin(), diss_rate.end(), diss_rate.begin(), detail::common__turbulence__tke<real_t>(lambda));
    }
  };
};
