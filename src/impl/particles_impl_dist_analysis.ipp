// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <boost/math/tools/minima.hpp>
#include <libcloudph++/common/detail/toms748.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    // init
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::dist_analysis_sd_conc(
      const common::unary_function<double> &n_of_lnrd_stp,
      const impl::n_t sd_conc,
      const real_t dt
    )
    {
      if(opts_init.rd_min >= 0 && opts_init.rd_max >= 0) // user-defined bin edges
      {
        real_t rd_min = opts_init.rd_min, rd_max = opts_init.rd_max;

        multiplier = log(rd_max / rd_min)
          / sd_conc
          * dt
          * (n_dims == 0
            ? dv[0]
            : (opts_init.dx * opts_init.dy * opts_init.dz)
          );

        log_rd_min = log(rd_min);
        log_rd_max = log(rd_max);
      }
      else if (opts_init.rd_min < 0 && opts_init.rd_max < 0) // automatic detection of bin edges
      {
        // probing the spectrum to find rd_min-rd_max range     
        // when analysing distro for source, multiplier takes into account that
        // the distribution is assumed to represent number of particles created per unit of time! 
        // TODO: document that

        // values to start the search 
        real_t rd_min = config.rd_min_init, rd_max = config.rd_max_init;

        bool found_optimal_range = false;
        while (!found_optimal_range)
        {
          multiplier = log(rd_max / rd_min)
            / sd_conc
            * dt
            * (n_dims == 0
              ? dv[0]
              : (opts_init.dx * opts_init.dy * opts_init.dz)
            );

          log_rd_min = log(rd_min);
          log_rd_max = log(rd_max);

          impl::n_t
            n_min = n_of_lnrd_stp(log_rd_min) * multiplier,
            n_max = n_of_lnrd_stp(log_rd_max) * multiplier;

          if (rd_min == config.rd_min_init && n_min != 0)
            throw std::runtime_error(detail::formatter() << "Initial dry radii distribution is non-zero (" << n_min << ") for rd_min_init (" << config.rd_min_init <<")");
          if (rd_max == config.rd_max_init && n_max != 0)
            throw std::runtime_error(detail::formatter() << "Initial dry radii distribution is non-zero (" << n_max << ") for rd_max_init (" << config.rd_max_init <<")");

          if      (n_min == 0) rd_min *= 1.01;
          else if (n_max == 0) rd_max /= 1.01;
          else found_optimal_range = true;
        }
      }
      else assert(false && "opts_init.rd_min * opts_init.rd_max < 0");
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::dist_analysis_const_multi(
      const common::unary_function<double> &n_of_lnrd_stp
    )
    {
      if(opts_init.rd_min >= 0 && opts_init.rd_max >= 0) // user-defined bin edges
      {
        log_rd_min = log(opts_init.rd_min);
        log_rd_max = log(opts_init.rd_max);
      }
      else if (opts_init.rd_min < 0 && opts_init.rd_max < 0) // automatic detection of bin edges
      {
        // TODO: add same sanity check as above
        //       how does brent algorithm work for functions with multiple minima??
        std::pair<real_t, real_t> init_distr_max; // [ln(position of distribution's maximum), -function value at maximum]
        boost::uintmax_t n_iter = config.n_iter;
        init_distr_max = boost::math::tools::brent_find_minima(detail::eval_and_mul<real_t>(n_of_lnrd_stp, -1), log(config.rd_min_init), log(config.rd_max_init), 200, n_iter); // bits = 200 to make algorithm choose max precision available

        real_t init_dist_bound_value = -init_distr_max.second / config.threshold; // value of the distribution at which we bind it
        n_iter = config.n_iter;
        // TODO: it could be written more clearly by creating an object detail::eval_and_oper<real_t>(*n_of_lnrd_stp, -init_dist_bound_value, 1), but for some reason it doesnt give the correct values
        log_rd_min = 
          common::detail::toms748_solve(
            detail::eval_and_add<real_t>(n_of_lnrd_stp, -init_dist_bound_value),
            real_t(log(config.rd_min_init)), init_distr_max.first,
            detail::eval_and_add<real_t>(n_of_lnrd_stp, -init_dist_bound_value)(real_t(log(config.rd_min_init))),
            detail::eval_and_add<real_t>(n_of_lnrd_stp, -init_dist_bound_value)(init_distr_max.first),
            config.eps_tolerance, n_iter
          );

        n_iter = config.n_iter;
        log_rd_max = 
          common::detail::toms748_solve(
            detail::eval_and_add<real_t>(n_of_lnrd_stp, -init_dist_bound_value),
            init_distr_max.first, real_t(log(config.rd_max_init)),
            detail::eval_and_add<real_t>(n_of_lnrd_stp, -init_dist_bound_value)(init_distr_max.first),
            detail::eval_and_add<real_t>(n_of_lnrd_stp, -init_dist_bound_value)(real_t(log(config.rd_max_init))),
            config.eps_tolerance, n_iter
          );
      }
      else assert(false && "opts_init.rd_min * opts_init.rd_max < 0");
    }
  };
};
