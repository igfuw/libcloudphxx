#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
    using common::unary_function;

    template <typename real_t>
    struct kappa_soluble_fraction_t {
          real_t kappa;
          real_t soluble_fraction; // volume fraction of the soluble part (0-1)

      kappa_soluble_fraction_t(real_t kappa_, real_t soluble_fraction_)
      : kappa(kappa_), soluble_fraction(soluble_fraction_)
      {}

      bool operator<(const kappa_soluble_fraction_t &other) const
      {
        if (kappa != other.kappa) return kappa < other.kappa;
        return soluble_fraction < other.soluble_fraction;
      }
        };

    // initial dry sizes of aerosol
    // defined with a distribution
    // uses shared_ptr to make opts_init copyable
    template<typename real_t>
    using dry_distros_t = std::map<
      kappa_soluble_fraction_t<real_t>,              // (kappa, soluble_fraction); dry_distros defines total dry radius
      std::shared_ptr<unary_function<real_t>> // n(ln(rd)) @ STP; alternatively it's n(ln(rd)) independent of rhod if aerosol_independent_of_rhod=true
    >;

    // defined with a size-number pair
    template<typename real_t>
    using dry_sizes_t = std::map<
      kappa_soluble_fraction_t<real_t>, // (kappa, soluble_fraction); dry_sizes defines total dry radius
      std::map<real_t,           // radius [m]
        std::pair<real_t, int>   // STP_concentration [1/m^3], number of SD that represent this radius kappa and concentration
      >
    >;

    // similar, but for sources of aerosols after initialization
    template<typename real_t>
    using src_dry_distros_t = std::map<
      kappa_soluble_fraction_t<real_t>,              // (kappa, soluble_fraction); src_dry_distros defines total dry radius
      std::tuple<std::shared_ptr<unary_function<real_t>>, int, int> // 1st: n(ln(rd)) @ STP created per second; alternatively it's n(ln(rd)) independent of rhod if aerosol_independent_of_rhod=true; 2nd: sd_conc for this distribution ; 3rd: supstp for this aerosol (interval in timesteps beween addition of these aerosols)
    >;

    // defined with a size-number pair
    template<typename real_t>
    using src_dry_sizes_t = std::map<
      kappa_soluble_fraction_t<real_t>, // (kappa, soluble_fraction); src_dry_sizes defines total dry radius
      std::map<real_t,           // radius [m]
        std::tuple<real_t, int, int>   // STP_concentration [1/m^3] created per second, number of SD that represent this radius kappa and concentration, supstp
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
