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
    void particles_t<real_t, device>::impl::init_hskpng_ncell()
    {
      T.resize(n_cell);
      p.resize(n_cell);
      RH.resize(n_cell); 
      eta.resize(n_cell); 

      count_ijk.resize(n_cell);
      count_num.resize(n_cell);
      count_mom.resize(n_cell);
      count_n = 0;
    }
  };
};
