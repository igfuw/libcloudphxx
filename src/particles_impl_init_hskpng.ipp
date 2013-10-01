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
    void particles_t<real_t, device>::impl::init_hskpng()
    {
      // memory allocation
      i.resize(n_part); // 
      j.resize(n_part); //  > TODO: are they needed at all?
      k.resize(n_part); // 
      ijk.resize(n_part);

      vt.resize(n_part);
      thrust::fill(vt.begin(), vt.end(), 0); // so that it may be safely used in condensation before first update

      T.resize(n_cell);
      p.resize(n_cell);
      r.resize(n_cell); 
      RH.resize(n_cell); 
      eta.resize(n_cell); 

      sorted_id.resize(n_part);
      sorted_ijk.resize(n_part);

      count_ijk.resize(n_cell);
      count_num.resize(n_cell);
      count_mom.resize(n_cell);
      count_n = 0;
    }
  };
};
