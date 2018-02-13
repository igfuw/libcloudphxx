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
    particles_t<real_t, device>::particles_t(opts_init_t<real_t> opts_init, int n_x_tot)
    {
      int rank, size;

      // handle MPI init
#if defined(USE_MPI)
      detail::mpi_init(MPI_THREAD_SINGLE, rank, size); 
#else
      rank = 0;
      size = 1;
      // throw an error if ran with mpi, but not compiled for mpi
      if (
        // mpich
        std::getenv("PMI_RANK") != NULL ||
        // openmpi
        std::getenv("OMPI_COMM_WORLD_RANK") != NULL ||
        // lam
        std::getenv("LAMRANK") != NULL ||
        // mvapich2
        std::getenv("MV2_COMM_WORLD_RANK") != NULL
      ) throw std::runtime_error("mpirun environment variable detected but libcloudphxx was compiled with MPI disabled");
#endif
      std::pair<detail::bcond_t, detail::bcond_t> bcond;
      if(size > 1)
        bcond = std::make_pair(detail::distmem_mpi, detail::distmem_mpi);
      else
        bcond = std::make_pair(detail::sharedmem, detail::sharedmem);

      // use the desired GPU card, TODO: remove it? can be done using CUDA_VISIBLE_DEVICES
#if defined(__NVCC__)
      if(opts_init.dev_id >= 0)
        cudaSetDevice(opts_init.dev_id);
#endif
      if(opts_init.dev_count < 2) // no distmem
        n_x_tot = opts_init.nx;

      // create impl instance
      pimpl.reset(new impl(opts_init, bcond, rank, size, n_x_tot));
      this->opts_init = &pimpl->opts_init;
      pimpl->sanity_checks();

      // init output map to 0
      for(int i=0; i < chem_all+2; ++i)
        pimpl->output_puddle[static_cast<output_t>(i)] = 0.;
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