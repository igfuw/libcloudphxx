// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <thrust/iterator/transform_iterator.h>
#include <libcloudph++/common/detail/toms748.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      // The condition for immersion freezing
      template<class real_t>
      class immersion_freeze_cond
      {
      public:
        BOOST_GPU_ENABLED
        bool operator()(const thrust::tuple<real_t, real_t, real_t> &tpl)  // tpl is a tuple of 3 elements: (T_freeze, ambient T, ambient RH)
        {
          if (thrust::get<0>(tpl) >=  thrust::get<1>(tpl) && thrust::get<2>(tpl) >= real_t(1))  // returns true if T_freeze >= ambient T and ambient RH >= 1
            return true;
          else
            return false;
        };
      };
    };

    // Immersion freezing
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::ice_nucl() {
      thrust::replace_if(ice.begin(), ice.end(),                        // Replacing values of ice with 1 if immersion_freeze_cond is satisfied.
        thrust::make_zip_iterator(
          thrust::make_tuple(                                           // Creating a zip iterator to access multiple vectors:
            T_freeze.begin(),                                               // freezing temperature for each droplet
            thrust::make_permutation_iterator(T.begin(), ijk.begin()),      // ambient temperature
            thrust::make_permutation_iterator(RH.begin(), ijk.begin())      // ambient RH
          )
        ),
        detail::immersion_freeze_cond<real_t>(),
        real_t(1)
      );
    }

    //TODO: latent heat of freezing
  }
}
