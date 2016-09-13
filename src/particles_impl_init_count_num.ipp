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

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_count_num_const_multi(
      const common::unary_function<real_t> *n_of_lnrd_stp, 
      const real_t &ratio)
    {
      const real_t precision = 1e-4; // size of bins in ln(radius) when calculating roots, integral, CDF
      const real_t integral = detail::integrate(*n_of_lnrd_stp, log_rd_min, log_rd_max, precision);

      // number of SDs per cell under STP conditions
      real_t multiplier = round(integral
        * ratio
        / real_t(opts_init.sd_const_multi)
        * opts_init.dx
        * opts_init.dy
        * opts_init.dz);

      namespace arg = thrust::placeholders;
      using common::earth::rho_stp;
      // initialize number of SDs in cells taking into account differences in rhod
      // and that not all Eulerian cells are fully covered by Lagrangian domain (round to int)
      thrust::transform(rhod.begin(), rhod.end(), dv.begin(), count_num.begin(),
        (multiplier * arg::_1 / real_t(rho_stp<real_t>() / si::kilograms * si::cubic_metres) * arg::_2 / (opts_init.dx * opts_init.dy * opts_init.dz) + real_t(0.5))
      );
    }
  };
};
