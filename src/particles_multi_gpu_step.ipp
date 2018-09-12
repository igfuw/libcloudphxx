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
    // time-stepping methods
    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::step_sync(
      const opts_t<real_t> &opts,
      arrinfo_t<real_t> th,
      arrinfo_t<real_t> rv,
      const arrinfo_t<real_t> rhod,
      const arrinfo_t<real_t> courant_1,
      const arrinfo_t<real_t> courant_2,
      const arrinfo_t<real_t> courant_3,
      std::map<enum chem_species_t, arrinfo_t<real_t> > ambient_chem
    )
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::step_sync, opts, th, rv, rhod, courant_1, courant_2, courant_3, ambient_chem);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::sync_in(
      arrinfo_t<real_t> th,
      arrinfo_t<real_t> rv,
      const arrinfo_t<real_t> rhod,
      const arrinfo_t<real_t> courant_1,
      const arrinfo_t<real_t> courant_2,
      const arrinfo_t<real_t> courant_3,
      std::map<enum chem_species_t, arrinfo_t<real_t> > ambient_chem
    )
    {
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::sync_in, th, rv, rhod, courant_1, courant_2, courant_3, ambient_chem);
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::step_cond(
      const opts_t<real_t> &opts,
      arrinfo_t<real_t> th,
      arrinfo_t<real_t> rv,
      std::map<enum chem_species_t, arrinfo_t<real_t> > ambient_chem
    )
    {
printf("mCUDA step_cond start\n");
      pimpl->mcuda_run(&particles_t<real_t, CUDA>::step_cond, opts, th, rv, ambient_chem);

//      detail::barrier_t barrier(this->opts_init->dev_count);
//printf("mCCUDA step_cond post barrier\n");
//if(this->opts_init->dev_id == 0)
//{
//printf("step cond dev id %d\n", this->opts_init->dev_id);
//      // do step async on each device
//      gpuErrchk(cudaSetDevice(this->opts_init->dev_id));
//      pimpl->particles[this->opts_init->dev_id]->step_cond(opts, th, rv, ambient_chem);
//printf("step cond dev id %d done\n", this->opts_init->dev_id);
//}
//        barrier.wait();
//if(this->opts_init->dev_id == 1)
//{
//printf("step cond dev id %d\n", this->opts_init->dev_id);
//      // do step async on each device
//      gpuErrchk(cudaSetDevice(this->opts_init->dev_id));
//      pimpl->particles[this->opts_init->dev_id]->step_cond(opts, th, rv, ambient_chem);
//printf("step cond dev id %d done\n", this->opts_init->dev_id);
//}
//        barrier.wait();
//if(this->opts_init->dev_id == 2)
//{
//printf("step cond dev id %d\n", this->opts_init->dev_id);
//      // do step async on each device
//      gpuErrchk(cudaSetDevice(this->opts_init->dev_id));
//      pimpl->particles[this->opts_init->dev_id]->step_cond(opts, th, rv, ambient_chem);
//printf("step cond dev id %d done\n", this->opts_init->dev_id);
//}
//        barrier.wait();
//if(this->opts_init->dev_id == 3)
//{
//printf("step cond dev id %d\n", this->opts_init->dev_id);
//      // do step async on each device
//      gpuErrchk(cudaSetDevice(this->opts_init->dev_id));
//      pimpl->particles[this->opts_init->dev_id]->step_cond(opts, th, rv, ambient_chem);
//printf("step cond dev id %d done\n", this->opts_init->dev_id);
//}
//        barrier.wait();
    }

    template <typename real_t>
    void particles_t<real_t, multi_CUDA>::step_async(
      const opts_t<real_t> &opts
    )
    {
      // sanity checks
      if(opts.rcyc)
        throw std::runtime_error("Particle recycling can't be used in the multi_CUDA backend (it would consume whole memory quickly");

      // cuda streams and events to control asynchronous copies
      // note: storing them in particles_multi_t caused errors
      // on program exit
      std::vector<cudaStream_t> streams(this->opts_init->dev_count);
      std::vector<cudaEvent_t> events(this->opts_init->dev_count);
      // if made a member of particles_t<multi_CUDA...> it causes errors on program exit...
      detail::barrier_t barrier(this->opts_init->dev_count);

      // run on many GPUs
      std::vector<std::thread> threads;
      for (int i = 0; i < this->opts_init->dev_count; ++i)
      {
        threads.emplace_back(
          &particles_t<real_t, multi_CUDA>::impl::step_async_and_copy, pimpl.get(), opts, i, std::ref(streams), std::ref(events), std::ref(barrier)
        );
      }
      for (auto &th : threads) th.join();
    }
  };
};
