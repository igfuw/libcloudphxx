#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
    using common::unary_function;

    // initial dry sizes of aerosol
    // defined with a distribution
    // uses shared_ptr to make opts_init copyable
    template<typename real_t>
    using dry_distros_t = std::map<
//      real_t,
      std::pair<real_t, real_t>,              // (kappa, ice)
      std::shared_ptr<unary_function<real_t>> // n(ln(rd)) @ STP; alternatively it's n(ln(rd)) independent of rhod if aerosol_independent_of_rhod=true
    >;

    // defined with a size-number pair
    template<typename real_t>
    using dry_sizes_t = std::map<
//      real_t,
      std::pair<real_t, real_t>, // (kappa, ice)
      std::map<real_t,           // radius [m]
        std::pair<real_t, int>   // STP_concentration [1/m^3], number of SD that represent this radius kappa and concentration
      >
    >;
  };
};
