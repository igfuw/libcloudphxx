#include <thrust/logical.h>

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    bool particles_t<real_t, device>::impl::perparticle_drw2_all_converged() { 
      const auto &unconverged_mask = sstp_cond_unconverged_mask_gp->get();
      return thrust::none_of(
        unconverged_mask.begin(), unconverged_mask.end(),
        cuda::std::identity()
      );
    }
  };
};