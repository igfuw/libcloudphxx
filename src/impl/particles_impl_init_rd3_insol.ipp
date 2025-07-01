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
    void particles_t<real_t, device>::impl::init_rd3_insol(
      const real_t &rd_insol // fixed rd of insoluble particles
    )
    {
      // half of super-particles contain insoluble substance
      thrust::fill(rd3_insol.begin() + n_part_old,
        rd3_insol.begin() + n_part_old + (rd3_insol.size() - n_part_old) / 2,
        rd_insol * rd_insol * rd_insol);

      // the other half contain soluble substances (rd_insol = 0)
      thrust::fill(rd3_insol.begin() + n_part_old + (rd3_insol.size() - n_part_old) / 2,
        rd3_insol.end(),
        real_t(0));
    }
  };
};

//TODO: 0.05 of insoluble particles should be inactive (for inactive the freezing temperature is -38)

