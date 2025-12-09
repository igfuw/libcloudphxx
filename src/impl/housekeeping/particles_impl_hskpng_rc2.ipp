// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief initialisation routine for super droplets
  */

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_approximate_rc2_invalid() // approximated for the temperature opts_init.rc2_T
    {
      if(sstp_cond_act == 1 || !allow_sstp_cond)
        return;

      namespace arg = thrust::placeholders;

      thrust::transform_if(
        rd3.begin(), rd3.end(), // input - 1st arg
        thrust::make_zip_iterator(make_tuple(
          kpa.begin(), 
          thrust::make_constant_iterator<real_t>(opts_init.rc2_T + 273.15)
        )),                                   // input - 2nd arg 
        rc2.begin(),                          // condition argument
        rc2.begin(),                          // output
        detail::rw3_cr<real_t>(),             // op
        arg::_1 == real_t(detail::invalid)    // condition
      );
    }
  };
};

