// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief initialisation routine for super droplets
  */

#include <libcloudph++/lgrngn/ice_nucleation.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_T_freeze() {

      using namespace libcloudphxx::lgrngn::ice_nucleation;

      // random numbers between [0,1] for sampling
      auto u01g = tmp_device_real_part.get_guard();
      thrust_device::vector<real_t> &u01 = u01g.get();
      rand_u01(u01, n_part_to_init);

      thrust::transform(
        thrust::make_zip_iterator(thrust::make_tuple(rd2_insol.begin() + n_part_old, u01.begin() + n_part_old)),
        thrust::make_zip_iterator(thrust::make_tuple(rd2_insol.end(),   u01.end())),
        T_freeze.begin() + n_part_old,
        T_freeze_CDF_inv_functor<real_t>(opts_init.inp_type)
      );
    }
  };
};

