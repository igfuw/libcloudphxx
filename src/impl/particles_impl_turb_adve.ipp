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
    void particles_t<real_t, device>::impl::turb_adve()
    {   
      thrust_device::vector<real_t> * vel_pos_a[] = {&x, &y, &z};
      std::vector<thrust_device::vector<real_t>*> vel_pos(&vel_pos_a[0], &vel_pos_a[0]+n_dims);

      thrust_device::vector<real_t> * vel_turbs_vctrs_a[] = {&up, &wp, &vp};
      std::vector<thrust_device::vector<real_t>*> vel_turbs_vctrs(&vel_turbs_vctrs_a[0], &vel_turbs_vctrs_a[0]+n_dims);

      namespace arg = thrust::placeholders;

      for(int i=0; i<n_dims; ++i)
      {
        thrust::transform(
          vel_pos[i]->begin(),
          vel_pos[i]->end(),
          vel_turbs_vctrs[i]->begin(),
          vel_pos[i]->begin(),
          arg::_1 + arg::_2 * opts_init.dt
        );
      }
    }
  };
};
