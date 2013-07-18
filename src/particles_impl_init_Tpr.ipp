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
    template <typename real_t, int device>
    void particles<real_t, device>::impl::init_Tpr()
    {
      // memory allocation
      rhod.resize(n_cell);
      rhod_th.resize(n_cell);
      rhod_rv.resize(n_cell);
      T.resize(n_cell);
      p.resize(n_cell);
      r.resize(n_cell); 
    }
  };
};
