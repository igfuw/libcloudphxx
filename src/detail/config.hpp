#pragma once
#include <libcloudph++/common/detail/toms748.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      // TODO: replace it with a namespace?
      //       this struct is instantiated twice
      template<class real_t>
      struct config
      {
         const uintmax_t n_iter = 100;      // number of iterations of the toms748 root-finding

        // precision of toms748
         const common::detail::eps_tolerance<real_t> eps_tolerance = 
          common::detail::eps_tolerance<real_t>(sizeof(real_t) * 8 / 4);

         const real_t bin_precision = 1e-4; // size of bins in ln(radius) when calculating roots, integral, CDF
         const real_t threshold = 1e4;      // in constant multiplicity mode, min and max initial radii are
                                           // defined by n_lnrd(r_min/max) = (maximum of n_lnrd) / threshold
         const real_t rd_min_init = 1e-11, 
                     rd_max_init = 1e-3;   // bounding values for the initial dry radius distro
      };
    };
  };
};
