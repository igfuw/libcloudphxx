// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  * @brief initialisation routine for super droplets
  */

namespace libcloudphxx
{
    namespace lgrngn
    {
        template <typename real_t, backend_t device>
        void particles_t<real_t, device>::impl::init_a_c_rho_ice()
        {
            // filling a and c with zeros initially
            thrust::fill(ice_a.begin() + n_part_old, ice_a.end(), real_t(0));
            thrust::fill(ice_c.begin() + n_part_old, ice_c.end(), real_t(0));
            thrust::fill(ice_rho.begin() + n_part_old, ice_rho.end(), real_t(0));
        }
    };
};