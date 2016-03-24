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
    void particles_t<real_t, device>::impl::unpack_n(const int &n_copied)
    {
      n_part_old = n_part;
      n_part += n_copied;
      n.resize(n_part);
      thrust::copy(in_n_bfr.begin(), in_n_bfr.begin() + n_copied, n.begin() + n_part_old);
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::unpack_real(const int &n_copied)
    {
      // real_t vectors copied
      thrust_device::vector<real_t> * real_t_vctrs[] = {&rd3, &rw2, &kpa, &vt, &x, &z, &y};
      const int real_vctrs_count =
        n_dims == 3 ? 7 :
          n_dims == 2 ? 6 : 5;

      for(int i = 0; i < real_vctrs_count; ++i)
      {
        real_t_vctrs[i]->resize(n_part);
        thrust::copy( in_real_bfr.begin() + i * n_copied, in_real_bfr.begin() + (i+1) * n_copied, real_t_vctrs[i]->begin() + n_part_old);
      }
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::flag_lft()
    {
      thrust::copy(
        thrust::make_constant_iterator<n_t>(0),
        thrust::make_constant_iterator<n_t>(0) + lft_count,
        thrust::make_permutation_iterator(n.begin(), lft_id.begin())
      );
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::flag_rgt()
    {
      thrust::copy(
        thrust::make_constant_iterator<n_t>(0),
        thrust::make_constant_iterator<n_t>(0) + rgt_count,
        thrust::make_permutation_iterator(n.begin(), rgt_id.begin())
      );
    }
  };
};

