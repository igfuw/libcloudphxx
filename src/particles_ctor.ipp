/// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief Thrust-based CPU/GPU particle-tracking logic for Lagrangian microphysics
  */

#include <thrust/host_vector.h>
#include <thrust/iterator/constant_iterator.h>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_algebra.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_operations.hpp>
#include <boost/numeric/odeint/external/thrust/thrust_resize.hpp>

#include <map>

namespace libcloudphxx
{
  namespace lgrngn
  {
    // ctor
    template <typename real_t, backend_t device>
    particles_t<real_t, device>::particles_t(const opts_init_t<real_t> &opts_init, const int &n_x_bfr, int n_x_tot) 
    {
#if defined(__NVCC__)
      if(opts_init.dev_id >= 0)
        cudaSetDevice(opts_init.dev_id);
#endif
      if(opts_init.dev_count < 2) // no distmem
        n_x_tot = opts_init.nx;

      pimpl.reset(new impl(opts_init, n_x_bfr, n_x_tot));

      this->opts_init = &pimpl->opts_init;
      pimpl->sanity_checks();

      // init output map to 0
      for(int i=0; i < chem_all+2; ++i)
        pimpl->output_puddle[static_cast<libcloudphxx::common::output_t>(i)] = 0.;
    }

    // dtor
    template <typename real_t, backend_t device>
    particles_t<real_t, device>::~particles_t() {};

    // outbuf
    template <typename real_t, backend_t device>
    real_t *particles_t<real_t, device>::outbuf() 
    {
      pimpl->fill_outbuf();
      // restore the count_num and count_ijk arrays
      pimpl->hskpng_count();
      return &(*(pimpl->tmp_host_real_cell.begin()));
    }
  };
};
