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
      template <typename real_t>
      struct remote
      {
        real_t lcl, rmt;

        remote(real_t lcl, real_t rmt) : lcl(lcl), rmt(rmt) {}

        BOOST_GPU_ENABLED
        real_t operator()(real_t x)
        {
          real_t res = rmt + x - lcl;
          if(res == rmt) res = nextafter(res, real_t(0.)); // in single precision, we used to get x=x1, TODO: does it call CUDA's nextafter when used by multi_CUDA?
          return res;
        }
      };
    };
 
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::pack_n_lft()
    {
      assert(out_n_bfr.size() >= lft_count);
      assert(in_n_bfr.size() >= lft_count);
      thrust::copy(
        thrust::make_permutation_iterator(n.begin(), lft_id.begin()),
        thrust::make_permutation_iterator(n.begin(), lft_id.begin()) + lft_count,
        out_n_bfr.begin()
      );
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::pack_n_rgt()
    {
      assert(out_n_bfr.size() >= rgt_count);
      assert(in_n_bfr.size() >= rgt_count);
      thrust::copy(
        thrust::make_permutation_iterator(n.begin(), rgt_id.begin()),
        thrust::make_permutation_iterator(n.begin(), rgt_id.begin()) + rgt_count,
        out_n_bfr.begin()
      );
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::pack_real_lft()
    {
      assert(out_real_bfr.size() >= lft_count * distmem_real_vctrs.size());
      assert(in_real_bfr.size() >= lft_count * distmem_real_vctrs.size());
      for(int i = 0; i < distmem_real_vctrs.size(); ++i)
        thrust::copy(
          thrust::make_permutation_iterator(distmem_real_vctrs[i]->begin(), lft_id.begin()),
          thrust::make_permutation_iterator(distmem_real_vctrs[i]->begin(), lft_id.begin()) + lft_count,
          out_real_bfr.begin() + i * lft_count
        );
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::pack_real_rgt()
    {
      assert(out_real_bfr.size() >= rgt_count * distmem_real_vctrs.size());
      assert(in_real_bfr.size() >= rgt_count * distmem_real_vctrs.size());
      for(int i = 0; i < distmem_real_vctrs.size(); ++i)
        thrust::copy(
          thrust::make_permutation_iterator(distmem_real_vctrs[i]->begin(), rgt_id.begin()),
          thrust::make_permutation_iterator(distmem_real_vctrs[i]->begin(), rgt_id.begin()) + rgt_count,
          out_real_bfr.begin() + i * rgt_count
        );
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::bcnd_remote_lft(const real_t &x0, const real_t &x1)
    {
      thrust::transform(
        thrust::make_permutation_iterator(x.begin(), lft_id.begin()),
        thrust::make_permutation_iterator(x.begin(), lft_id.begin()) + lft_count,
        thrust::make_permutation_iterator(x.begin(), lft_id.begin()), // in place
        detail::remote<real_t>(x0, x1)
      );
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::bcnd_remote_rgt(const real_t &x1, const real_t &x0)
    {
      thrust::transform(
        thrust::make_permutation_iterator(x.begin(), rgt_id.begin()),
        thrust::make_permutation_iterator(x.begin(), rgt_id.begin()) + rgt_count,
        thrust::make_permutation_iterator(x.begin(), rgt_id.begin()), // in place
        detail::remote<real_t>(x1, x0)
      );
    }
  };
};

