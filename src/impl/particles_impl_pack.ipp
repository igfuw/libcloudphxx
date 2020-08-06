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

      auto it = distmem_real_vctrs.begin();

      while (it != distmem_real_vctrs.end())
      {
        thrust::copy(
          thrust::make_permutation_iterator((*it)->begin(), lft_id.begin()),
          thrust::make_permutation_iterator((*it)->begin(), lft_id.begin()) + lft_count,
          out_real_bfr.begin() + std::distance(distmem_real_vctrs.begin(), it) * lft_count
        );
        it++;
      }
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::pack_real_rgt()
    {
      assert(out_real_bfr.size() >= rgt_count * distmem_real_vctrs.size());
      assert(in_real_bfr.size() >= rgt_count * distmem_real_vctrs.size());

      auto it = distmem_real_vctrs.begin();

      while (it != distmem_real_vctrs.end())
      {
        thrust::copy(
          thrust::make_permutation_iterator((*it)->begin(), rgt_id.begin()),
          thrust::make_permutation_iterator((*it)->begin(), rgt_id.begin()) + rgt_count,
          out_real_bfr.begin() + std::distance(distmem_real_vctrs.begin(), it) * rgt_count
        );
        it++;
      }
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
#if !defined(NDEBUG)
      auto min_it = thrust::min_element(
        thrust::make_permutation_iterator(x.begin(), lft_id.begin()),
        thrust::make_permutation_iterator(x.begin(), lft_id.begin()) + lft_count);
      if(*min_it < x0)
      {
        std::cerr << "x (" << *min_it << ")  < x0 (" << x0 << ") after adjustment for distmem copy, potentially SD moved by more than one process/GPU domain size"
        assert(0);
      }
#endif
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
#if !defined(NDEBUG)
      auto max_it = thrust::max_element(
        thrust::make_permutation_iterator(x.begin(), lft_id.begin()),
        thrust::make_permutation_iterator(x.begin(), lft_id.begin()) + lft_count);
      if(*max_it >= x1)
      {
        std::cerr << "x (" << *max_it << ")  >= x1 (" << x1 << ") after adjustment for distmem copy, potentially SD moved by more than one process/GPU domain size"
        assert(0);
      }
#endif
    }
  };
};

