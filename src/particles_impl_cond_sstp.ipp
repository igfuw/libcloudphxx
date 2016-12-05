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
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::cond_sstp(
      const real_t &dt,
      const real_t &RH_max
    ) {   
      // prerequisite
      hskpng_sort(); 
      // particle's local change in rv
      thrust_device::vector<real_t> &pdrv(tmp_device_real_part4);
      // -rw3_old
      thrust::transform(
        thrust::make_transform_iterator(rw2.begin(), detail::rw2torw3<real_t>()),
        thrust::make_transform_iterator(rw2.end(), detail::rw2torw3<real_t>()),
        pdrv.begin(),
        thrust::negate<real_t>()
      );

      // vector for each particle's T
      thrust_device::vector<real_t> &Tp(tmp_device_real_part3);

      // calc Tp
      thrust::transform(
        sstp_tmp_th.begin(), sstp_tmp_th.end(), // input - first arg
        sstp_tmp_rh.begin(),                    // input - second arg
        Tp.begin(),                             // output
        detail::common__theta_dry__T<real_t>() 
      );  

      // calculating drop growth in a timestep using backward Euler 
      thrust::transform(
        rw2.begin(), rw2.end(),         // input - 1st arg (zip not as 1st arg not to write zip.end()
        thrust::make_zip_iterator(      // input - 2nd arg
          thrust::make_tuple(
            sstp_tmp_rh.begin(),
            sstp_tmp_rv.begin(),
            Tp.begin(),
            // particle-specific p
            thrust::make_transform_iterator(
              thrust::make_zip_iterator(
                thrust::make_tuple(
                  sstp_tmp_rh.begin(),
                  sstp_tmp_rv.begin(),
                  Tp.begin()
              )),
              detail::common__theta_dry__p<real_t>()
            ),
            // particle-specific RH
            thrust::make_transform_iterator(
              thrust::make_zip_iterator(
                thrust::make_tuple(
                  sstp_tmp_rh.begin(),
                  sstp_tmp_rv.begin(),
                  Tp.begin()
              )),
              detail::RH<real_t>()
            ),
            // particle-specific eta
            thrust::make_transform_iterator(
              Tp.begin(),
              detail::common__vterm__visc<real_t>()
            ),
            rd3.begin(),
            kpa.begin(),
            vt.begin()
          )
        ), 
        rw2.begin(),                    // output
        detail::advance_rw2<real_t>(dt, RH_max, &config)
      );

      // calc rw3_new - rw3_old
      thrust::transform(
        thrust::make_transform_iterator(rw2.begin(), detail::rw2torw3<real_t>()),
        thrust::make_transform_iterator(rw2.end(), detail::rw2torw3<real_t>()),
        pdrv.begin(),
        pdrv.begin(),
        thrust::plus<real_t>()
      );

      // calc - 4/3 * pi * rho_w * n * (rw3_new - rw3_old) / (dV * rhod)
      thrust::transform(
        pdrv.begin(), pdrv.end(),                  // input - 1st arg
        thrust::make_zip_iterator(thrust::make_tuple(
          sstp_tmp_rh.begin(),                                        // rhod
          n.begin(),                                                  // n
          thrust::make_permutation_iterator(dv.begin(), ijk.begin()) // dv
        )),
        pdrv.begin(),                             // output
        detail::rw3diff2drv<real_t>(
          - common::moist_air::rho_w<real_t>() / si::kilograms * si::cubic_metres
          * real_t(4./3) * pi<real_t>(), n_dims
        )
      );  

      // apply change in rv to sstp_tmp_rv
      update_pstate(sstp_tmp_rv, pdrv);

      // calc particle-specific change in th based on pdrv
      thrust::transform(
        thrust::make_zip_iterator(thrust::make_tuple(  
          pdrv.begin(),       //  
          Tp.begin(),         // dth = drv * d_th_d_rv(T, th)
          sstp_tmp_th.begin() //  
        )), 
        thrust::make_zip_iterator(thrust::make_tuple(  
          pdrv.end(),       //  
          Tp.end(),         // dth = drv * d_th_d_rv(T, th)
          sstp_tmp_th.end() //  
        )), 
        pdrv.begin(), // in-place
        detail::dth<real_t>()
      );

      // apply change in th to sstp_tmp_th
      update_pstate(sstp_tmp_th, pdrv);
    }
  };  
};
