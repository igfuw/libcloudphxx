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
      return (bcond.first == detail::distmem_mpi || bcond.second == detail::distmem_mpi);
    }
    template <typename real_t, backend_t device>
    bool particles_t<real_t, device>::impl::distmem_cuda(
    )
    {
      return (detail::bcond_is_distmem_cuda(bcond.first) || detail::bcond_is_distmem_cuda(bcond.second));
    }
    template <typename real_t, backend_t device>
    bool particles_t<real_t, device>::impl::distmem(
    )
    {
      return (distmem_mpi() || distmem_cuda());
    }
  };
};
