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
    // init
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::dist_analysis(
      const common::unary_function<real_t> *n_of_lnrd_stp, // TODO: kappa-spectrum map
      const impl::n_t sd_conc,
      const real_t dt
    )
    {
      // probing the spectrum to find rd_min-rd_max range     
      // when analysing distro for source, multiplier takes into account that
      // the distribution is assumed to represent number of particles created per unit of time! 
      // TODO: document that

      // values to start the search 
      const real_t rd_min_init = 1e-11, rd_max_init = 1e-3;
      real_t rd_min = rd_min_init, rd_max = rd_max_init;

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
          n_min = (*n_of_lnrd_stp)(log_rd_min) * multiplier,
          n_max = (*n_of_lnrd_stp)(log_rd_max) * multiplier;

        if (rd_min == rd_min_init && n_min != 0)
          throw std::runtime_error("Initial dry radii distribution is non-zero for rd_min_init");
        if (rd_max == rd_max_init && n_max != 0)
          throw std::runtime_error("Initial dry radii distribution is non-zero for rd_max_init");

        if      (n_min == 0) rd_min *= 1.1;
        else if (n_max == 0) rd_max /= 1.1;
        else found_optimal_range = true;
      }
    };
  };
};
