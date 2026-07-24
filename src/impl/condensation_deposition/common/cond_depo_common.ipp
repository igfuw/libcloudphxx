// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

// #include <thrust/iterator/transform_iterator.h>
#include <libcloudph++/common/maxwell-mason.hpp>
#include <libcloudph++/common/kappa_koehler.hpp>
#include <libcloudph++/common/kelvin_term.hpp>
#include <libcloudph++/common/transition_regime.hpp>
#include <libcloudph++/common/ventil.hpp>
#include <libcloudph++/common/mean_free_path.hpp>
#include <libcloudph++/common/detail/toms748.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template<class real_t>
      struct massdiff2drv 
      {
        real_t mlt;
        int n_dims;
        
        BOOST_GPU_ENABLED
        massdiff2drv(const real_t &mlt, const int &n_dims):
          mlt(mlt), n_dims(n_dims) {}
 
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &massdiff, const thrust::tuple<real_t, real_t, real_t> &tpl) noexcept
        {
          if(n_dims > 0)
            return mlt * massdiff * thrust::get<1>(tpl) / thrust::get<0>(tpl) / thrust::get<2>(tpl);
          else // for parcel setup use 1/rhod instead of dv, dv will be updated in hskpng_Tpr in async
            return mlt * massdiff * thrust::get<1>(tpl);
        }
      };

       template<class real_t, int power>
       struct rw2torwX
       {
         BOOST_GPU_ENABLED
         real_t operator()(const real_t &rw2) noexcept
         {
 #if !defined(__NVCC__)
           using std::pow;
 #endif
           return pow(rw2, real_t(power) / real_t(2));
         }
       };

       template<class real_t>
       struct rw2torwX<real_t, 3>
       {
         BOOST_GPU_ENABLED
         real_t operator()(const real_t &rw2) noexcept
         {
 #if !defined(__NVCC__)
           using std::sqrt;
 #endif
           return rw2 * sqrt(rw2);
         }
       };

       template<class real_t>
       struct rw2torwX<real_t, 2>
       {
         BOOST_GPU_ENABLED
         real_t operator()(const real_t &rw2) noexcept
         {
           return rw2;
         }
       };        
    };
  };  
};
