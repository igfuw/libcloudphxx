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
      template <typename real_t>
      struct chem_vol_fun
      { // calculate drop volume
        const real_t pi;

        // ctor (pi() is not a __device__ function...)
        chem_vol_fun() :
          pi(boost::math::constants::pi<real_t>())
        {}

        BOOST_GPU_ENABLED
        real_t operator()(const real_t &rw2) const
        {
#if !defined(__NVCC__)
	  using std::pow;
#endif
          return real_t(4./3) * pi * (pow(rw2, real_t(3./2)));
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::chem_vol_ante()
    {   
      if (opts_init.chem_switch == false) throw std::runtime_error("all chemistry was switched off");

      //calculate new drop volumes (to be used in Henrys law)
      thrust_device::vector<real_t> &V(tmp_device_real_part);

      thrust::transform(
        rw2.begin(), rw2.end(),         // input
        V.begin(),                      // output 
        detail::chem_vol_fun<real_t>()  // op
      );
    }

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::chem_vol_post()
    {   
      if (opts_init.chem_switch == false) throw std::runtime_error("all chemistry was switched off");

      //save the current drop volume in V_old (to be used in the next step for Henrys law)
      thrust_device::vector<real_t> &V(tmp_device_real_part);
      thrust_device::vector<real_t> &V_old(tmp_device_real_part_V_old);
    
      thrust::copy(V.begin(), V.end(), V_old.begin());
    }
  };  
};
