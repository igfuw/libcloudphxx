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
    template <typename real_t, int device>
    void particles<real_t, device>::sync_e2l(
      typename parent_t::arrinfo_t rhod_th,
      typename parent_t::arrinfo_t rhod_rv,
      typename parent_t::arrinfo_t rhod     // defaults to {NULL, NULL}
    )  
    {
std::cerr << "\n\n SYNC_E2L \n\n";
    }

    template <typename real_t, int device>
    void particles<real_t, device>::sync_l2e(
      typename parent_t::arrinfo_t rhod_th,
      typename parent_t::arrinfo_t rhod_rv
    )  
    {
std::cerr << "\n\n SYNC_L2E \n\n";
    }
  };
};
