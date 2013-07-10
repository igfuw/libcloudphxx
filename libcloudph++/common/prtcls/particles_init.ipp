// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <iostream>

#include "particles.hpp"

#include "detail/urand.hpp"
#include "detail/thrust.hpp"

namespace libcloudphxx
{
  namespace common
  {
    namespace prtcls
    {
      // init
      template <typename real_t, int device>
      void particles<real_t, device>::init(
        real_t rd_min, real_t rd_max
      )
      {
        using namespace thrust::placeholders;

std::cerr << "init called, calling urand..." << std::endl;
        // tossing random numbers [0,1]
        pimpl->urand(pimpl->n_part);

        // shifting from [0,1] to [log(rd_min),log(rd_max)] and storing into rd3
        // pimpl->rd3 temporarily means logarithm of radius!
        thrust::transform(
          pimpl->u01.begin(), 
          pimpl->u01.end(), 
          pimpl->rd3.begin(), 
          log(rd_min) + _1 * (log(rd_max) - log(rd_min)) // TODO: does not work!
        );
 
        //thrust::transform();
debug::print(pimpl->rd3);

std::cerr << "done." << std::endl;
      }
    };
  };
};
