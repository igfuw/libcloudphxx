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
    template <typename real_t>
    particles_t<real_t, multi_CUDA>::particles_t(opts_init_t<real_t> _opts_init) 
    {
      pimpl.reset(new impl(_opts_init));
  
      // make opts_init point to global opts init
      this->opts_init = &(pimpl->glob_opts_init);
    }

    // dtor
    template <typename real_t>
    particles_t<real_t, multi_CUDA>::~particles_t() {}

    // TODO: what about MPI with other backends?
    // initialisation 
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::init(
      const arrinfo_t<real_t> th,
      const arrinfo_t<real_t> rv,
      const arrinfo_t<real_t> rhod,
      const arrinfo_t<real_t> courant_1,
      const arrinfo_t<real_t> courant_2,
      const arrinfo_t<real_t> courant_3,
      const std::map<enum chem_species_t, const arrinfo_t<real_t> > ambient_chem
    )
    {
      pimpl->mcuda_run(
        &particles_t<real_t, CUDA>::init,
        th, rv, rhod, courant_1, courant_2, courant_3, ambient_chem
      );
      // exchange domains using mpi; has to be done sequentially here as MPI isn't good with calls from many threads per node; UPDATE: actually, MPI implementations used by libmpdata support such calls
      // TODO: hence, move this exchange to some &particles_t<real_t, multi_CUDA>::init ?
      pimpl->mcuda_run(
        &particles_t<real_t, CUDA>::xchng_domains
      );
/*
      if(this->opts_init->dev_count > 1)
      {
        // first node receives first 
        if(pimpl->particles[0]->pimpl->mpi_rank==0)
        {
          pimpl->particles[0]->pimpl->xchng_domains();
          pimpl->particles[this->opts_init->dev_count-1]->pimpl->xchng_domains();
        }
        else  // other nodes send first
        {
          pimpl->particles[this->opts_init->dev_count-1]->pimpl->xchng_domains();
          pimpl->particles[0]->pimpl->xchng_domains();
        }
      }
      else
        pimpl->particles[0]->pimpl->xchng_domains();
*/
    }
  };
};
