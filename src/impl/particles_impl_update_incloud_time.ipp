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
    namespace detail
    {
      template <class real_t>
      struct incloud_time_updater
      {
        const real_t dt;

        incloud_time_updater(const real_t &dt):
          dt(dt) {}

        template<class tup_t>
        BOOST_GPU_ENABLED
        void operator()(tup_t tup)
        {
          if(thrust::get<0>(tup) > thrust::get<1>(tup)) // rw > rc
            thrust::get<2>(tup) += dt;
          else
            thrust::get<2>(tup) = 0;

        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::update_incloud_time(const real_t &dt)
    {   
      // tmp vector to store crit radius of each SD
      thrust_device::vector<real_t> &rc2(tmp_device_real_part);

      // computing rc2
      thrust::transform(
        rd3.begin(), rd3.end(), // input - 1st arg
        thrust::make_zip_iterator(make_tuple(
          kpa.begin(), 
          thrust::make_permutation_iterator(
            T_ref.begin(),
            ijk.begin_ref()
          )
        )),                                   // input - 2nd arg 
        rc2.begin(),                          // output
        detail::rw3_cr<real_t>()              // op
      );

      // increase time if activated, set to 0 if not
      thrust::for_each_n(
        thrust::make_zip_iterator(make_tuple(
          rw2.begin(),
          rc2.begin(),
          incloud_time.begin()
        )),
        n_part,
        detail::incloud_time_updater<real_t>(dt)
      );
    }
  };  
};
