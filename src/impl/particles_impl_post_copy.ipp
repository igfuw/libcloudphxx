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
    void particles_t<real_t, device>::impl::post_copy(const opts_t<real_t> &opts)
    {
      // recycling out-of-domain/invalidated particles 
      if(opts.rcyc)
        rcyc();
      // if we do not recycle, we should remove them
      else
        hskpng_remove_n0();  

      // updating particle->cell look-up table
      hskpng_ijk();

      // updating count_ijk and count_num
      hskpng_count();
    }
  };
};
