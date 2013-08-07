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
    // init
    template <typename real_t, int device>
    void particles<real_t, device>::impl::init_wet()
    {
      // memory allocation
      rw2.resize(n_part);
      kpa.resize(n_part);
    }
  };
};
