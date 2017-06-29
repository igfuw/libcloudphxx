// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <thrust/sequence.h>
#include <algorithm>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace{
    void comb(int N, int K, thrust::host_vector<int> &vec, int off)
    {
        std::string bitmask(K, 1); // K leading 1's
        bitmask.resize(N, 0); // N-K trailing 0's
    
        // print integers and permute bitmask
        do {
            for (int i = 0; i < N; ++i) // [0..N-1] integers
            {
                if (bitmask[i]) vec.push_back(i+off);
            }
        } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
    }
    };

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
        rand_un(n_part);

        // sorting the sequence with the random key
        thrust::sort_by_key(
          un.begin(), un.end(),
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

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_shuffle_allperm_and_sort()
    {   
      hskpng_count(); // get number of SDs per cell

      // helper host vector with count_nums
      thrust::host_vector<int> h_count_num(n_cell);
      thrust::copy(count_num.begin(), count_num.end(), h_count_num.begin());

      namespace arg = thrust::placeholders;
      // count number of SD pairs
      thrust::transform( 
        count_num.begin(),
        count_num.end(),
        count_num.begin(),
        arg::_1 * (arg::_1 - 1) / 2
      );
      // get total number of pairs
      n_tot_col_pairs = thrust::reduce(count_num.begin(), count_num.end(), 0);

      // helper host vector...
      thrust::host_vector<int> perm;//(2*n_tot_col_pairs);
      thrust::host_vector<int> h_n_sd_b4(n_cell);
      thrust::exclusive_scan(h_count_num.begin(), h_count_num.end(), h_n_sd_b4.begin());

      for(int i=0; i < n_cell; ++i) 
      {
        comb(h_count_num[i], 2, perm, h_n_sd_b4[i]);
      }
//debug::print(perm);
      thrust_device::vector<int> d_perm(2*n_tot_col_pairs);
      thrust::copy(perm.begin(), perm.end(), d_perm.begin());

      // get sorted id's
      hskpng_sort_helper(false); 
//debug::print(sorted_id);
      col_pairs.resize(2*n_tot_col_pairs);
      thrust::copy(
        thrust::make_permutation_iterator(sorted_id.begin(), d_perm.begin()),
        thrust::make_permutation_iterator(sorted_id.begin(), d_perm.begin()) + 2*n_tot_col_pairs,
        col_pairs.begin()
      );
      sorted_ijk_col.resize(2*n_tot_col_pairs);
      thrust::copy(
        thrust::make_permutation_iterator(ijk.begin(), col_pairs.begin()),
        thrust::make_permutation_iterator(ijk.begin(), col_pairs.begin()) + 2*n_tot_col_pairs,
        sorted_ijk_col.begin()
      );
//debug::print(col_pairs);





      // resize the sorted_id, that will store indices of colliding SDs
      //sorted_id.resize(2*n_tot_col_pairs);
    }
  };  
};
