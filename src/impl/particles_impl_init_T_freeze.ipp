// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief initialisation routine for super droplets
  */

#include "../../include/libcloudph++/common/ice_nucleation.hpp"

namespace libcloudphxx
{
  namespace lgrngn
  {
    template <typename real_t, backend_t device>
  void particles_t<real_t, device>::impl::init_T_freeze(
    const int rng_seed
    )
    {
      using namespace libcloudphxx::common::ice_nucleation;

      const INP_t inp_type = INP_t::mineral; //TODO: INP type as argument, to support different partcle types

      thrust::transform(
        rd3_insol.begin() + n_part_old,
        rd3_insol.end(),
        T_freeze.begin() + n_part_old,
        T_freeze_functor<real_t>(inp_type, rng_seed)
      );
    }
  };
};

