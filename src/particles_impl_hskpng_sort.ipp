// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <thrust/sequence.h>

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, int device>
    void particles<real_t, device>::impl::hskpng_sort()
    {   
      // e.g. after shuffling
      if (sorted) return;

      // filling-in sorted_id with a sequence
      thrust::sequence(sorted_id.begin(), sorted_id.end());

      // making a copy of ijk
      thrust::copy(
        ijk.begin(), ijk.end(), // from
        sorted_ijk.begin()      // to
      );

      // sorting sorted_ijk and sorted_id
      thrust::sort_by_key(
	sorted_ijk.begin(), sorted_ijk.end(), // keys
	sorted_id.begin()                     // values
      );

      // flagging that particles are now sorted
      sorted = true;
    }   
  };  
};
