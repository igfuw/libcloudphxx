// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief initialisation routine for super droplets
  */

#include <libcloudph++/common/ice_nucleation.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::init_T_freeze() {

      using namespace libcloudphxx::common::ice_nucleation;

      const INP_t inp_type = INP_t::mineral; //TODO: INP type as argument, to support different partcle types

      // random numbers between [0,1] for sampling
      rand_u01(n_part_to_init);

      thrust::transform(
        thrust::make_zip_iterator(thrust::make_tuple(rd2_insol.begin() + n_part_old, u01.begin() + n_part_old)),
        thrust::make_zip_iterator(thrust::make_tuple(rd2_insol.end(),   u01.end())),
        T_freeze.begin() + n_part_old,
        T_freeze_CDF_inv_functor<real_t>(inp_type)
      );
    }
  };
};

