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
    // calc the SGS turbulent velocity component
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_sgs_vel_ST_periodic(const real_t &dt, const bool only_vertical)
    {   
      if(only_vertical) throw std::runtime_error("ST_periodic only_vertical doesn't work yet.");

      namespace arg = thrust::placeholders;

      // progress time in the model
      ST.update_time(dt);

//    zero-out SGS velocities
      thrust::fill(up.begin(), up.end(), real_t(0));
      thrust::fill(vp.begin(), vp.end(), real_t(0));
      thrust::fill(wp.begin(), wp.end(), real_t(0));

      // calc new SGS velocities
      auto zip_pos_vel = 
        thrust::make_zip_iterator(thrust::make_tuple(
          x.begin(),
          y.begin(),
          z.begin(),
          up.begin(),
          vp.begin(),
          wp.begin()
      ));

      ST.calc_vel(zip_pos_vel, n_part);
    };
  };
};
