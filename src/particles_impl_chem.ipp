// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <libcloudph++/common/detail/bisect.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t>
      struct chem_minfun
      {
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &m_H) const
        {
          return 0;
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::chem(const real_t &dt)
    {   
std::cerr << "@particles_t::impl::chem()" << std::endl;
      // 1/3: equilibrium stuff

      real_t tol = 44;
      real_t range = common::detail::bisect(
        detail::chem_minfun<real_t>(),
        real_t(22), // min -> pure water
        real_t(44), // max -> TODO
        tol
      ); 

      // 2/3: non-equilibrium stuff
      // TODO


      // 3/3: recomputing dry radii
      {
        using namespace thrust::placeholders;
        // TODO: using namespace for S_VI
        thrust::transform(
          che[S_VI].begin(), che[S_VI].end(),                      // input
          rd3.begin(),                                             // output
          (real_t(3./4) / pi<real_t>() / opts_init.chem_rho) * _1  // op
        );
      };
    }
  };  
};
