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
      struct sedi_with_vt_and_subsidence
      {
        const real_t dt;

        sedi_with_vt_and_subsidence(real_t _dt): dt(_dt) {}

        template<class tuple_t>
        BOOST_GPU_ENABLED
        real_t operator()(real_t z, tuple_t tpl) const
        {
          return z - dt * (thrust::get<0>(tpl) + thrust::get<1>(tpl));
        };
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::sedi()
    {   
      namespace arg = thrust::placeholders;
 
      // settling due to sedimentation + large-scale subsidence
      thrust::transform(
        z.begin(), z.end(),                    // position
        thrust::make_zip_iterator(
          thrust::make_tuple(
            vt.begin(),                                                    // terminal velocity 
            thrust::make_permutation_iterator(w_LS.begin(), k.begin())     // large-scale subsidence velocity
          )
        ),
        z.begin(),                         // output
        detail::sedi_with_vt_and_subsidence<real_t>(opts_init.dt)
      );
    }
  };  
};
