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
    namespace detail
    {
      template <typename real_t, typename n_t>
      struct scale_factor
      {
        BOOST_GPU_ENABLED
        real_t operator()(const n_t &n)
        {
          // see section 5.1.3 in Shima et al. 2009
          return real_t((n*(n-1))/2) / (n/2); 
        }
      };
    };

    template <typename real_t, int device>
    void particles<real_t, device>::impl::coal(const real_t &dt)
    {   
      // prerequisites
      hskpng_shuffle_and_sort(); // to get random neighbours by default
      hskpng_count();            // no. of particles per cell 
      
      // placing scale_factors in count_mom (of size count_n!)
      thrust::transform(
        count_num.begin(), count_num.begin() + count_n, // input - 1st arg
        count_mom.begin(),                              // output
        detail::scale_factor<real_t, n_t>()
      );

      // laying out scale factor onto ijk grid
      thrust_device::vector<real_t> &scl(tmp_device_real_cell);
      thrust::fill(scl.begin(), scl.end(), real_t(0));
      thrust::copy(
        count_mom.begin(),                    // input - begin
        count_mom.begin() + count_n,          // input - end
        thrust::make_permutation_iterator(    // output
          scl.begin(),                        // data
          count_ijk.begin()                   // permutation
        )
      );  

      // colliding
      
    }
  };  
};
