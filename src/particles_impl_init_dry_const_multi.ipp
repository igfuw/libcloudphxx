// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief initialisation routine for super droplets
  */

#include <iostream>
#include <algorithm>

#include "detail/thrust.hpp"
#include "detail/functors_host.hpp"

#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust/extrema.h>
#include <thrust/binary_search.h>

#include <boost/math/tools/minima.hpp>

#include <libcloudph++/common/earth.hpp>

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

      // calculate cumulative distribution function
      template<typename real_t, typename vec_t>
      void calc_CDF(const common::unary_function<real_t> &fun, const real_t &min, const real_t &max, const real_t &bin_size, vec_t &vec)
      {
        const int n = (max - min) / bin_size + 1; //no of points at which cdf will be calculated
        vec.resize(n);

        namespace arg = thrust::placeholders;
        // fill vec with fun values at n points
        thrust::transform(
          thrust::make_transform_iterator(thrust::make_counting_iterator(0), min + bin_size * arg::_1),
          thrust::make_transform_iterator(thrust::make_counting_iterator(n), min + bin_size * arg::_1),
          vec.begin(), eval_and_oper<real_t>(fun, 1));

        // calculate CDF
        thrust::inclusive_scan(vec.begin(), vec.end(), vec.begin());
       
        // normalize CDF     
        thrust::transform(vec.begin(), vec.end(), vec.begin(), arg::_1 / vec.back());
      }
    };

    // init
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_dry_const_multi(
      const real_t kappa,
      const common::unary_function<real_t> *n_of_lnrd_stp // TODO: kappa-spectrum map
    )
    {
      // values to start the search 
      const real_t rd_min_init = 1e-11, rd_max_init = 1e-3;
      real_t rd_min = rd_min_init, rd_max = rd_max_init;

      boost::uintmax_t max_iter = 1e6; // max number of iterations when searching for maxima/roots
      const real_t precision = 1e-4; // size of bins in ln(radius) when calculating roots, integral, CDF

      // analysing initial distribution
      // has to have only single maximum, no local minima
      std::pair<real_t, real_t> init_distr_max; // [ln(position of distribution's maximum), -function value at maximum]
      init_distr_max = boost::math::tools::brent_find_minima(detail::eval_and_oper<real_t>(*n_of_lnrd_stp, -1), log(rd_min), log(rd_max), 32, max_iter);
  
      real_t init_dist_bound_value = -init_distr_max.second / 1e4; // value of the distribution at which we bind it

      rd_min = 
        common::detail::bisect<libcloudphxx::lgrngn::detail::eval_and_oper<real_t>, real_t>(
          detail::eval_and_oper<real_t>(*n_of_lnrd_stp, -init_dist_bound_value, 1), 
          real_t(log(rd_min_init)), init_distr_max.first, 
          precision
        );

      rd_max = 
        common::detail::bisect<libcloudphxx::lgrngn::detail::eval_and_oper<real_t>, real_t>(
          detail::eval_and_oper<real_t>(*n_of_lnrd_stp, -init_dist_bound_value, 1), 
          init_distr_max.first, real_t(log(rd_max_init)),
          precision
        );

      const real_t integral = detail::integrate(*n_of_lnrd_stp, rd_min, rd_max, precision);

      // calculate cumulative distribution function
      thrust::host_vector<real_t> cdf;

      detail::calc_CDF(*n_of_lnrd_stp, rd_min, rd_max, precision, cdf);

      // number of SDs per cell under STP conditions
      real_t multiplier = round(integral 
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

      n_part = thrust::reduce(count_num.begin(), count_num.end());

      // memory allocation
      rd3.resize(n_part);
      n.resize(n_part);
      kpa.resize(n_part); 
      tmp_device_real_part.resize(n_part);
      tmp_device_n_part.resize(n_part);

      // filling kappas
      thrust::fill(kpa.begin(), kpa.end(), kappa);
 
      // filling multiplicities
      thrust::fill(n.begin(), n.end(), opts_init.sd_const_multi);

      // tossing random numbers [0,1] for dry radii
      rand_u01(n_part);

      // rd3 temporarily means logarithm of radius!
      thrust_device::vector<real_t> &lnrd(rd3);

      thrust::host_vector<real_t> host_u01(n_part); 
      thrust::copy(u01.begin(), u01.end(), host_u01.begin());
      thrust::host_vector<real_t> host_lnrd(n_part); 
      
      // sample ln(rd) from the distribution with the inverse transform sampling method
      thrust::upper_bound(cdf.begin(), cdf.end(), host_u01.begin(), host_u01.end(), host_lnrd.begin());
      thrust::copy(host_lnrd.begin(), host_lnrd.end(), lnrd.begin());
      thrust::transform(lnrd.begin(), lnrd.end(), lnrd.begin(), rd_min + arg::_1 * precision);

      // converting rd back from logarithms to rd3
      thrust::transform(
        lnrd.begin(),
        lnrd.end(),
        rd3.begin(),
        detail::exp3x<real_t>()
      );
    }
  };
};
