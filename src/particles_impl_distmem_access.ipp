// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    bool particles_t<real_t, device>::impl::distmem_mpi(
    )
    {
      return bcond.first == detail::distmem_mpi || bcond.second == detail::distmem_mpi;
    }
    template <typename real_t, backend_t device>
    bool particles_t<real_t, device>::impl::distmem_cuda(
    )
    {
      return bcond.first == detail::distmem_cuda || bcond.second == detail::distmem_cuda;
    }
    template <typename real_t, backend_t device>
    bool particles_t<real_t, device>::impl::distmem(
    )
    {
      return distmem_mpi() || distmem_cuda();
    }
  };
};
