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

      if (!shuffle)
      {
	// making a copy of ijk
	thrust::copy(
	  ijk.begin(), ijk.end(), // from
	  sorted_ijk.begin()      // to
	);
      }
      else
      {
        // generating a random sorting key
        rand_u01(n_part);

        // sorting the sequence with the random key
        thrust::sort_by_key(
          u01.begin(), u01.end(),
          sorted_id.begin()
        );

        // permuting sorted_ijk accordingly
        thrust::copy(
          thrust::make_permutation_iterator(ijk.begin(), sorted_id.begin()), // input - begin
          thrust::make_permutation_iterator(ijk.end(),   sorted_id.end()  ), // input - end
          sorted_ijk.begin()                                                 // output
        );
      }

      // sorting sorted_ijk and sorted_id
      thrust::sort_by_key(
	sorted_ijk.begin(), sorted_ijk.end(), // keys
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
