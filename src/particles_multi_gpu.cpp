// contains definitions of members of particles_t specialized for multiple GPUs

namespace libcloudphxx
{
  namespace lgrngn
  {

    // constructor
    template <typename real_t>
    particles_t<real_t, multi_CUDA>::particles_t(const opts_init_t<real_t> &)
    {
//    int dev_count; // number of GPUs used
//    std::vector<particles_t<real_t, CUDA> *> particles; // pointer to particles_t on each GPU
    }

    // initialisation 
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::init(
      const arrinfo_t<real_t> th,
      const arrinfo_t<real_t> rv,
      const arrinfo_t<real_t> rhod,
      const arrinfo_t<real_t> courant_1,
      const arrinfo_t<real_t> courant_2,
      const arrinfo_t<real_t> courant_3
    )
    {
    }

    // time-stepping methods
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::step_sync(
      const opts_t<real_t> &,
      arrinfo_t<real_t> th,
      arrinfo_t<real_t> rv,
      const arrinfo_t<real_t> courant_1,
      const arrinfo_t<real_t> courant_2,
      const arrinfo_t<real_t> courant_3,
      const arrinfo_t<real_t> rhod
    )
    {
    }

    template <typename real_t>
    real_t particles_t<real_t, multi_CUDA>::step_async(
      const opts_t<real_t> &
    )
    {
    }

    // diagnostic methods
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_sd_conc()
    {
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_dry_rng(
      const real_t &r_mi, const real_t &r_mx
    )
    {
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_wet_rng(
      const real_t &r_mi, const real_t &r_mx
    )
    {
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_dry_mom(const int &k)
    {
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_wet_mom(const int &k)
    {
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_wet_mass_dens(const real_t&, const real_t&)
    {
    }

    // ...
//</listing>

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_chem(const enum chem_species_t&)
    {
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_rw_ge_rc()
    {
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_RH_ge_Sc()
    {
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::diag_all()
    {
    }

  };
};
