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
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &diss_rate, const real_t &lambda)
        {
          assert(lambda > 0);
          return common::GA17_turbulence::tke(
            diss_rate * si::metres * si::metres / si::seconds / si::seconds / si::seconds,
            lambda * si::metres) / si::metres / si::metres * si::seconds * si::seconds;
        }
      };
    };
    // calc the SGS TKE, in place of dissipation rate
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_tke()
    {   
      namespace arg = thrust::placeholders;

      thrust::transform(
        diss_rate.begin(), diss_rate.end(),
        thrust::make_permutation_iterator(SGS_mix_len.begin(),   // profile of the SGS mixing length
          thrust::make_transform_iterator(                       // calculate vertical index from cell index
            thrust::make_counting_iterator<thrust_size_t>(0),
            arg::_1 % opts_init.nz
          )
        ),
        diss_rate.begin(),
        detail::common__turbulence__tke<real_t>());
    }
  };
};
