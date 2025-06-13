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
      template<class real_t>
      class immersion_freeze_cond
      {
        BOOST_GPU_ENABLED
        bool operator()(const auto &tpl)
        {
          if (thrust::get<0>(tpl) >=  thrust::get<1>(tpl) && thrust::get<2>(tpl) >= real_t(1))
            return true;
          else
            return false;
        };
      };
    };

  template<class real_t>
  void particles_t<real_t, device>::impl::nucleation() {
    thrust::replace_if(ice.begin(), ice.end(),
      thrust::make_zip_iterator(
        thrust::make_tuple(
          T_freeze.begin(),
          thrust::make_permutation_iterator(T.begin(), ijk.begin()),
          thrust::make_permutation_iterator(RH.begin(), ijk.begin())
        )
      ),
      detail::immersion_freeze_cond<real_t>(),
      real_t(1)
    );
  }
};
