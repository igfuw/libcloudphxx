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
    // init_xyz
    template <typename real_t, int device>
    void particles<real_t, device>::init_xyz()
    {
      //if (nx != 0) x.resize(n_part);
      //if (ny != 0) y.resize(n_part);
      //if (nz != 0) z.resize(n_part);

/*
      using namespace thrust::placeholders;

      // tossing random numbers [0,1] for dry radii
      pimpl->urand(pimpl->n_part);

	// shifting from [0,1] to [log(rd_min),log(rd_max)] and storing into rd3
	thrust::transform(
	  pimpl->u01.begin(), 
	  pimpl->u01.end(), 
	  lnrd.begin(), 
	  log(rd_min) + _1 * (log(rd_max) - log(rd_min)) 
	);
 
*/
    }
  };
};
