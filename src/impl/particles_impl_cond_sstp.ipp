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
      if(!const_p) // variable pressure
      {
        thrust::transform(
          sstp_tmp_th.begin(), sstp_tmp_th.end(), // input - first arg
          sstp_tmp_rh.begin(),                    // input - second arg
          Tp.begin(),                             // output
          detail::common__theta_dry__T_rhod<real_t>() 
        );  
      }
      else // external pressure profile
      {
        // T = dry2std(th_d, rv) * exner(p_tot)
        thrust::transform(
          sstp_tmp_th.begin(), sstp_tmp_th.end(),                      // input - first arg
          thrust::make_zip_iterator(thrust::make_tuple(
            sstp_tmp_rv.begin(),                                       // input - second arg 
            sstp_tmp_p.begin(),                                        // input - third arg
          )),
          Tp.begin(),                                                  // output
          detail::common__theta_dry__T_p<real_t>() 
        );
      }


      // calculating drop growth in a timestep using backward Euler 
      // TODO: these two function calls only differ in the way pressure is calculated,
      //       find a way to reuse it (e.g. std::bind?)
      // TODO2: in const_p we don't substep pressure
      if(!const_p)
      {
        // particle-specific pressure iterator, used twice
        // TODO: store the value somewhere?
        auto pressure_iter = 
          thrust::make_transform_iterator(
            thrust::make_zip_iterator(
              thrust::make_tuple(
                sstp_tmp_rh.begin(),
                sstp_tmp_rv.begin(),
                Tp.begin()
            )),
            detail::common__theta_dry__p<real_t>()
          );

        thrust::transform(
          rw2.begin(), rw2.end(),         // input - 1st arg (zip not as 1st arg not to write zip.end()
          thrust::make_zip_iterator(      // input - 2nd arg
            thrust::make_tuple(
              sstp_tmp_rh.begin(),
              sstp_tmp_rv.begin(),
              Tp.begin(),
              // particle-specific p
              pressure_iter,
              // particle-specific RH
              thrust::make_transform_iterator(
                thrust::make_zip_iterator(
                  thrust::make_tuple(
                    pressure_iter,
                    sstp_tmp_rv.begin(),
                    Tp.begin()
                )),
                detail::RH<real_t>(opts_init.RH_formula)
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
          detail::advance_rw2<real_t>(dt, RH_max)
        );
      }
      else
      {
        thrust::transform(
          rw2.begin(), rw2.end(),         // input - 1st arg (zip not as 1st arg not to write zip.end()
          thrust::make_zip_iterator(      // input - 2nd arg
            thrust::make_tuple(
              sstp_tmp_rh.begin(),
              sstp_tmp_rv.begin(),
              Tp.begin(),
              // particle-specific p
              sstp_tmp_p.begin(),
              // particle-specific RH
              thrust::make_transform_iterator(
                thrust::make_zip_iterator(
                  thrust::make_tuple(
                    sstp_tmp_p.begin(),
                    sstp_tmp_rv.begin(),
                    Tp.begin()
                )),
                detail::RH<real_t>(opts_init.RH_formula)
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
          detail::advance_rw2<real_t>(dt, RH_max)
        );
      }

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
