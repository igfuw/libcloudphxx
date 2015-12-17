// vim:filetype=cpp

/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief timestepping routine for super droplets
  */

namespace libcloudphxx
{
  namespace lgrngn
  {
    // some stuff to be done at the end of the step.
    // if using more than 1 GPU
    // has to be done after copy 
    template <typename real_t, backend_t device>
    real_t particles_t<real_t, device>::impl::step_finalize()
    {
      // recycling out-of-domain/invalidated particles 
      // (doing it here and not in async reduces the need for a second sort before diagnostics,
      // but also unneccesarily holds dyncore execution for a bit longer)
      thrust_size_t n_rcyc = 0;//pimpl->rcyc();
      // TODO: ! if we do not recycle, we should remove them to care for out-od-domain advection after sedimentation...

      // remove SDs with n = 0
      // for more than 2 GPUs, remove will be called from multi_gpu.cpp
      if(opts.sedi || opts.adve || opts.coal) 
        hskpng_remove_n0();  

      // updating particle->cell look-up table
      // (before advection and sedimentation so that their order does not matter,
      if (opts.adve || opts.sedi || n_rcyc)
        hskpng_ijk();
    }
  };
};
