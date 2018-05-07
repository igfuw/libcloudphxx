// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

// contains definitions of members of particles_t specialized for multiple GPUs
namespace libcloudphxx
{
  namespace lgrngn
  {
    // constructor
    template <typename real_t>
    particles_t<real_t, multi_CUDA>::particles_t(const opts_init_t<real_t> &_opts_init) 
    {
      pimpl.reset(new impl(_opts_init));
  
      // make opts_init point to global opts init
      this->opts_init = &(pimpl->glob_opts_init);
    }

    // dtor
    template <typename real_t>
    particles_t<real_t, multi_CUDA>::~particles_t() {}

    // initialisation 
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::init(
      const arrinfo_t<real_t> th,
      const arrinfo_t<real_t> rv,
      const arrinfo_t<real_t> rhod,
      const arrinfo_t<real_t> p,
      const arrinfo_t<real_t> p_d,
      const arrinfo_t<real_t> courant_1,
      const arrinfo_t<real_t> courant_2,
      const arrinfo_t<real_t> courant_3,
      const std::map<enum chem_species_t, const arrinfo_t<real_t> > ambient_chem
    )
    {
      pimpl->mcuda_run(
        &particles_t<real_t, CUDA>::init,
        th, rv, rhod, p, p_d, courant_1, courant_2, courant_3, ambient_chem
      );
    }
  };
};
