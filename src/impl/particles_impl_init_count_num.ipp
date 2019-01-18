// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */
#include <thrust/sequence.h>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      // calculate numerical integral with the trapezoidal rule
      // TODO: use thrust
      template<typename real_t>
      real_t integrate(const common::unary_function<real_t> &fun, const real_t &min, const real_t &max, const real_t &bin_size)
      {
        const int n = (max - min) / bin_size; //no of bins
        real_t integral = (fun(min) + fun(max)) / 2.;

        for(int i=1; i<n; ++i)
          integral += fun(min + i * bin_size);

        return integral * bin_size;
      }
    };

    // init number of SDs to be initialized per cell
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_count_num_sd_conc(const real_t &ratio)
    {
      thrust::fill(count_num.begin(), count_num.end(), ratio * opts_init.sd_conc);
    }

    // calculate number of droplets in a cell from concentration [1/m^3], taking into account cell volume and air density
    template <typename real_t, backend_t device>
    template <class arr_t>
    void particles_t<real_t, device>::impl::conc_to_number(arr_t &arr) 
    {
      assert(arr.size() == dv.size());

      namespace arg = thrust::placeholders;
      using common::earth::rho_stp;

      // cell volume
      thrust::transform(
        dv.begin(), dv.end(), 
        arr.begin(), 
        arr.begin(),
        arg::_2 * arg::_1 
      );

      // correct for density with respect to STP
      if(!opts_init.aerosol_independent_of_rhod)
        thrust::transform(
          rhod.begin(), rhod.end(), 
          arr.begin(), 
          arr.begin(),
          arg::_1 / real_t(rho_stp<real_t>() / si::kilograms * si::cubic_metres) * arg::_2  
        );
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_count_num_hlpr(const real_t &conc, const thrust_size_t &const_multi)
    {
      thrust_device::vector<real_t> &concentration(tmp_device_real_cell);
      thrust::fill(concentration.begin(), concentration.end(), conc);
      conc_to_number(concentration);

      namespace arg = thrust::placeholders;
      thrust::transform(
        concentration.begin(), concentration.end(),
        count_num.begin(),
        arg::_1 / const_multi + real_t(0.5)
      );
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_count_num_dry_sizes(const std::pair<real_t, int> &conc_multi)
    {
      thrust::fill(count_num.begin(), count_num.end(), conc_multi.second);
      //init_count_num_hlpr(conc_multi.first, conc_multi.second);
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_count_num_const_multi(
      const common::unary_function<real_t> &n_of_lnrd_stp
    )
    {
      const real_t integral = detail::integrate(n_of_lnrd_stp, log_rd_min, log_rd_max, config.bin_precision);
      init_count_num_hlpr(integral, opts_init.sd_const_multi);
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_count_num_const_multi(
      const common::unary_function<real_t> &n_of_lnrd_stp,
      const thrust_size_t &const_multi
    )
    {
      const real_t integral = detail::integrate(n_of_lnrd_stp, log_rd_min, log_rd_max, config.bin_precision);
      init_count_num_hlpr(integral, const_multi);
    }
  };
};
