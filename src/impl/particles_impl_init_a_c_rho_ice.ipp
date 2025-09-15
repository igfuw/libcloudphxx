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
            thrust::fill(a_ice.begin() + n_part_old, a_ice.end(), real_t(0));
            thrust::fill(c_ice.begin() + n_part_old, c_ice.end(), real_t(0));
            thrust::fill(rho_i.begin() + n_part_old, rho_i.end(), real_t(0));
        }
    };
};