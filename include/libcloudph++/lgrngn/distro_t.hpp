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
      std::tuple<real_t, real_t, unsigned long long, unsigned long long>, // kappa, soluble_fraction, sd_conc, sd_const_multi
      std::shared_ptr<unary_function<real_t>> // n(ln(rd)) @ STP; alternatively it's n(ln(rd)) independent of rhod if aerosol_independent_of_rhod=true
    >;

    // defined with a size-number pair
    template<typename real_t>
    using dry_sizes_t = std::map<
      std::tuple<real_t, real_t>, // (kappa, soluble_fraction); dry_sizes defines total dry radius
      std::map<real_t,           // radius [m]
        std::pair<real_t, int>   // STP_concentration [1/m^3], number of SD that represent this radius kappa and concentration
      >
    >;

    // similar, but for sources of aerosols after initialization
    template<typename real_t>
    using src_dry_distros_t = std::map<
      std::tuple<real_t, real_t, int, int, int>, // kappa, soluble_fraction, sd_conc, sd_const_multi, supstp
      std::shared_ptr<unary_function<real_t>> // n(ln(rd)) @ STP created per second; alternatively it's n(ln(rd)) independent of rhod if aerosol_independent_of_rhod=true
    >;

    // defined with a size-number pair
    template<typename real_t>
    using src_dry_sizes_t = std::map<
      std::tuple<real_t, real_t, int>, // kappa, soluble_fraction, supstp
      std::map<real_t,           // radius [m]
        std::pair<real_t, int>   // STP_concentration [1/m^3] created per second, number of SD that represent this radius kappa and concentration
      >
    >;

    // uses shared_ptr to make opts_init copyable
    // TODO: allow partial solubility (soluble_fraction < 1) for rlx_dry_distros (currently, its assumed to be 1)
    template<typename real_t>
    using rlx_dry_distros_t = std::unordered_map<
      real_t,                // kappa
      std::tuple<
        std::shared_ptr<unary_function<real_t>>, // n(ln(rd)) @ STP; alternatively it's n(ln(rd)) independent of rhod if aerosol_independent_of_rhod=true
        std::pair<real_t, real_t>, // kappa range of CCN considered to belong to this distribution, ranges of different members of the map need to be exclusive (TODO: add a check of this)
        std::pair<real_t, real_t>  // range of altitudes at which this relaxation acts
      >
    > ;
  };
};
