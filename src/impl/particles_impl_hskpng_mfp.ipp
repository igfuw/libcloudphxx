// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */


// Calculate mean free path used to calculate molecular correction for condensation
// NOTE: results are stored in tmp arrays

#include <libcloudph++/common/mean_free_path.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t>
      struct common__mean_free_path__lambda_D 
      {
        BOOST_GPU_ENABLED 
        real_t operator()(const real_t &T)
        {   
          return common::mean_free_path::lambda_D<real_t>(T  * si::kelvins) / si::meters;
        }   
      }; 

      template <typename real_t>
      struct common__mean_free_path__lambda_K 
      {
        BOOST_GPU_ENABLED 
        real_t operator()(const real_t &T, const real_t &p)
        {   
          return common::mean_free_path::lambda_K<real_t>(T  * si::kelvins, p * si::pascals) / si::meters;
        }   
      }; 
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_mfp()
    {   
      // tmp_device_real_cell cant be used here, it is used in cond.ipp
      thrust_device::vector<real_t> &lambda_D(opts_init.n_ref > 1 ? tmp_device_real_cell_ref  : tmp_device_real_cell1); 
      thrust_device::vector<real_t> &lambda_K(opts_init.n_ref > 1 ? tmp_device_real_cell_ref1 : tmp_device_real_cell2); 

      thrust::transform(T.begin(), T.end(), lambda_D.begin(), detail::common__mean_free_path__lambda_D<real_t>());
      thrust::transform(T.begin(), T.end(), p.begin(), lambda_K.begin(), detail::common__mean_free_path__lambda_K<real_t>());
    }
  };  
};
