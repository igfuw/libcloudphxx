// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#include <libcloudph++/common/GA17_turbulence.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template<class real_t>
      struct common__GA17__turbulence__tau
      {
        BOOST_GPU_ENABLED
        real_t operator()(const real_t &tke, const real_t &lambda)
        {
          assert(lambda > 0);
          return common::GA17_turbulence::tau(
            tke * si::metres * si::metres / si::seconds / si::seconds,
            lambda * si::metres) / si::seconds;
        }
      };

      template<class real_t>
      struct common__GA17__turbulence__update_sgs_vel
      {
        const quantity<si::time, real_t> dt;
        common__GA17__turbulence__update_sgs_vel(const real_t &dt):
          dt(dt * si::seconds){}

        template<class tpl_t>
        BOOST_GPU_ENABLED
        real_t operator()(tpl_t tpl)
        {
          return common::GA17_turbulence::update_sgs_vel(
            thrust::get<0>(tpl) * si::metres / si::seconds,
            thrust::get<1>(tpl) * si::seconds,
            dt,
            thrust::get<2>(tpl) * si::metres * si::metres / si::seconds / si::seconds,
            quantity<si::dimensionless, real_t>(thrust::get<3>(tpl))
          ) / si::metres * si::seconds;
        }
      };
    };
    // calc the SGS turbulent velocity component
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::hskpng_sgs_vel_GA17(const real_t &dt, const bool only_vertical)
    {   
      namespace arg = thrust::placeholders;

      thrust_device::vector<real_t> &tau(tmp_device_real_cell);
      thrust_device::vector<real_t> &tke(diss_rate); // should be called after hskpng_tke, which replaces diss_rate with tke

      thrust::transform(tke.begin(), tke.end(),
        thrust::make_permutation_iterator(SGS_mix_len.begin(),   // profile of the SGS mixing length
          thrust::make_transform_iterator(                       // calculate vertical index from cell index
            thrust::make_counting_iterator<thrust_size_t>(0),
            arg::_1 % opts_init.nz
          )
        ),
        tau.begin(),
        detail::common__GA17__turbulence__tau<real_t>()
      );

      thrust_device::vector<real_t> &r_normal(tmp_device_real_part);
      thrust_device::vector<real_t> * vel_turbs_vctrs_a[] = {&up, &wp, &vp};
      for(int i = (only_vertical ? 1 : 0); i < (only_vertical ? 2 : n_dims); ++i)
      {
        rng.generate_normal_n(r_normal, n_part); // generate a random number for wach particle with a normal distribution with mean 0 and std dev 1
        thrust::transform(
          thrust::make_zip_iterator(thrust::make_tuple(
            vel_turbs_vctrs_a[i]->begin(),
            thrust::make_permutation_iterator(tau.begin(), ijk.begin()),
            thrust::make_permutation_iterator(tke.begin(), ijk.begin()),
            r_normal.begin()
          )),
          thrust::make_zip_iterator(thrust::make_tuple(
            vel_turbs_vctrs_a[i]->begin(),
            thrust::make_permutation_iterator(tau.begin(), ijk.begin()),
            thrust::make_permutation_iterator(tke.begin(), ijk.begin()),
            r_normal.begin()
          )) + n_part,
          vel_turbs_vctrs_a[i]->begin(),
          detail::common__GA17__turbulence__update_sgs_vel<real_t>(dt)
        );
      }
    }
  };
};
