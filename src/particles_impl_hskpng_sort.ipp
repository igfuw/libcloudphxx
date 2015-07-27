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
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_sort_helper(bool shuffle)
    {   
      // filling-in sorted_id with a sequence
      thrust::sequence(sorted_id.begin(), sorted_id.end());

      // putting active SDs in front
      thrust::sort_by_key(
        sd_stat.begin(), sd_stat.end(),
        sorted_id.begin(),
        detail::active_first()
      );

      if (shuffle)
      {
        // generating a random sorting key
        rand_un(n_part);

        // sorting the active ids with the random key
        thrust::sort_by_key(
          un.begin(), un.begin() + n_part,
          sorted_id.begin()
        );
      }

      // permuting sorted_ijk accordingly
      thrust::copy(
        thrust::make_permutation_iterator(ijk.begin(), sorted_id.begin()), // input - begin
        thrust::make_permutation_iterator(ijk.begin(), sorted_id.begin() + n_part  ), // input - end
        sorted_ijk.begin()                                                 // output
      );

      // sorting sorted_ijk and sorted_id
      thrust::sort_by_key(
	sorted_ijk.begin(), sorted_ijk.begin() + n_part, // keys
	sorted_id.begin()                     // values
      );

      // flagging that particles are now sorted
      sorted = true;
    }   

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_sort()
    {   
      if (sorted) return; // e.g. after shuffling
      hskpng_sort_helper(false);
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_shuffle_and_sort()
    {   
      hskpng_sort_helper(true);
    }
  };  
};
