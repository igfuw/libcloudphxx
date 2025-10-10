/// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief Thrust-based CPU/GPU particle-tracking logic for Lagrangian microphysics
  */

#include <thrust/host_vector.h>
#include <thrust/iterator/constant_iterator.h>

#include <map>

namespace libcloudphxx
{
  namespace lgrngn
  {
    using namespace common::chem;

    // ctor
    template <typename real_t, backend_t device>
    particles_t<real_t, device>::particles_t(opts_init_t<real_t> opts_init, int n_x_tot)
    {
      int rank, size;

      // handle MPI init
#if defined(USE_MPI)
      detail::mpi_init(MPI_THREAD_FUNNELED, rank, size); 
#else
      rank = 0;
      size = 1;
      // throw an error if ran with mpi, but not compiled for mpi
      if ( ran_with_mpi() )
        throw std::runtime_error("libcloudph++: mpirun environment variable detected but libcloudphxx was compiled with MPI disabled");
#endif
      std::pair<detail::bcond_t, detail::bcond_t> bcond;
      if(size > 1)
      {
        if(!opts_init.open_side_walls) // periodic bcond in x
          bcond = std::make_pair(detail::distmem_mpi, detail::distmem_mpi);
        else // open bcond in x
        {
          if(rank == 0)
            bcond = std::make_pair(detail::open, detail::distmem_mpi);
          else if(rank == size-1)
            bcond = std::make_pair(detail::distmem_mpi, detail::open);
          else
            bcond = std::make_pair(detail::distmem_mpi, detail::distmem_mpi);
        }
      }
      else // only one process
      {
        if(!opts_init.open_side_walls) // periodic bcond in x
          bcond = std::make_pair(detail::sharedmem, detail::sharedmem);
        else
          bcond = std::make_pair(detail::open, detail::open);
      }

      // use the desired GPU card, TODO: remove it? can be done using CUDA_VISIBLE_DEVICES
#if defined(__NVCC__)
      if(opts_init.dev_id >= 0)
        gpuErrchk(cudaSetDevice(opts_init.dev_id));
#endif
      if(opts_init.dev_count < 2) // no distmem
        n_x_tot = opts_init.nx;

      // create impl instance
      pimpl.reset(new impl(opts_init, bcond, rank, size, n_x_tot));
      this->opts_init = &pimpl->opts_init;
      pimpl->sanity_checks();

      // init output map to 0
      for(int i=0; i < common::output_names.size(); ++i)
        pimpl->output_puddle[static_cast<common::output_t>(i)] = 0.;
    }

    // dtor
    template <typename real_t, backend_t device>
    particles_t<real_t, device>::~particles_t() {};

    // outbuf
    template <typename real_t, backend_t device>
    real_t *particles_t<real_t, device>::outbuf() 
    {
      auto outbuf_g = pimpl->tmp_host_real_cell.get_guard();
      thrust::host_vector<real_t> &outbuf = outbuf_g.get();

      pimpl->fill_outbuf(outbuf);
      // restore the count_num and count_ijk arrays
      pimpl->hskpng_count();
      return &(*(outbuf.begin()));
    }

    template <typename real_t, backend_t device>
    std::vector<real_t> particles_t<real_t, device>::get_attr(const std::string &attr_name) 
    {
      return std::move(pimpl->fill_attr_outbuf(attr_name));
    }
  };
};
