#include <libcloudph++/common/vterm.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template<class real_t>
      struct bin_mid
      {
        real_t ln_r_min, dlnr;
        bin_mid(const real_t &r_min, const real_t &r_max, const int &n_bin):
          ln_r_min(log(r_min)),
          dlnr(log(r_max / r_min) / n_bin)
          {}

        BOOST_GPU_ENABLED
        real_t operator()(const int &it)
        {
          return exp(ln_r_min + (it+0.5)*dlnr);
        }
      };

      struct vt_0
      {
        template<class real_t>
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &r)
        {
          return common::vterm::vt_beard77_v0(r * si::metres) / si::metres_per_second;
        }
      };
    };
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_vterm()
    {
      if(opts_init.terminal_velocity != vt_t::beard77fast) return;

      const int n_bin=4444;
      const real_t r_min = 1e-9;
      const real_t r_max = 3e-3; // 6mm is the max diameter Beard 1977 discusses
      vt_0.resize(n_bin);
      
      // calc mid r of each bin
      thrust::transform
      (
        thrust::make_counting_iterator<int>(0),
        thrust::make_counting_iterator<int>(0) + n_bin,
        vt_0.begin(),
        detail::bin_mid<real_t>(r_min, r_max, n_bin)
      );
      // calc vt_0
      thrust::transform
      (
        vt_0.begin(),
        vt_0.end(),
        vt_0.begin(),
        detail::vt_0()
      );
    }
  }
}
