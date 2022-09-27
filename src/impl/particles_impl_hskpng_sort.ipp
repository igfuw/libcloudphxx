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
//        thrust_device::vector<thrust_size_t> *ijk[] = {&ijk.get(), &ijk.get_ref()};
        thrust::copy(
          ijk.begin(), ijk.end(), // from
          sorted_ijk.begin()      // to
        );
        thrust::copy(
          ijk.begin_ref(), ijk.end_ref(), // from
          sorted_ijk.begin_ref()      // to
        );
      }
      else
      {
        // generating a random sorting key
        rand_un(n_part);

        // sorting the sequence with the random key
        thrust::sort_by_key(
          un.begin(), un.end(),
          sorted_id.begin()
        );

        // permuting sorted_ijk ref accordingly
        thrust::copy(
          thrust::make_permutation_iterator(ijk.begin_ref(), sorted_id.begin()), // input - begin
          thrust::make_permutation_iterator(ijk.end_ref(),   sorted_id.end()  ), // input - end
          sorted_ijk.begin_ref()                                                 // output
        );
      }

      // sorting sorted_ijk and sorted_id
      thrust::sort_by_key(
        sorted_ijk.begin_ref(), sorted_ijk.end_ref(), // keys
        sorted_id.begin()                             // values
      );

      // set sroted_ijk for normal grid
      thrust::copy(
        thrust::make_permutation_iterator(ijk.begin(), sorted_id.begin()),
        thrust::make_permutation_iterator(ijk.begin(), sorted_id.begin()) + n_part,
        sorted_ijk.begin()
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
