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
    void particles<real_t, device>::diag()
    {
std::cerr << "\n\n DIAG \n\n";
      // super-droplet concentration per grid cell
      pimpl->hskpng_count();
      
      thrust::fill(pimpl->tmp_host_real_cell.begin(), pimpl->tmp_host_real_cell.end(), 0);
      thrust::copy(
        pimpl->count_num.begin(),                  // input - begin
        pimpl->count_num.begin() + pimpl->count_n, // input - end
        thrust::make_permutation_iterator(         // output
          pimpl->tmp_host_real_cell.begin(),     
          pimpl->count_ijk.begin()
        )
      );
    }
  };
};
