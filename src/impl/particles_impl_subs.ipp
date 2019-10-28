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
    namespace detail
    {
      template<class real_t>
      struct subsidence
      {
        const real_t dt;

        subsidence(real_t _dt): dt(_dt) {}

        BOOST_GPU_ENABLED
        real_t operator()(real_t z, real_t v) const
        {
          return z - dt * v;
        };
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::subs()
    {   
      namespace arg = thrust::placeholders;
 
      // settling due to sedimentation + large-scale subsidence
      thrust::transform(
        z.begin(), z.end(),                    // position
        thrust::make_permutation_iterator(w_LS.begin(), k.begin()),     // large-scale subsidence velocity
        z.begin(),                         // output
        detail::subsidence<real_t>(opts_init.dt)
      );
    }
  };  
};
