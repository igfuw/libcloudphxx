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
    void particles_t<real_t, device>::impl::init_hskpng_ncell()
    {
      // memory allocation, p already resized in init_sync()
      T.resize(n_cell);
      RH.resize(n_cell);
      eta.resize(n_cell);
      if (opts_init.ice_switch)
        RH_i.resize(n_cell);

      count_ijk.resize(n_cell);
      count_num.resize(n_cell);
      count_mom.resize(n_cell);
      count_n = 0;

      // initialising device temporary arrays
      tmp_device_real_cell.resize(n_cell);
      tmp_device_size_cell.resize(n_cell);
      tmp_host_size_cell.resize(n_cell);
      tmp_host_real_cell.resize(n_cell);
      if(allow_sstp_cond && !opts_init.exact_sstp_cond)
      {
        sstp_tmp_rv.resize(n_cell);
        sstp_tmp_th.resize(n_cell);
        sstp_tmp_rh.resize(n_cell);
      }
      revp20.resize(n_cell, real_t(0));
      revp25.resize(n_cell, real_t(0));
      revp32.resize(n_cell, real_t(0));
      accr20.resize(n_cell, real_t(0));
      accr25.resize(n_cell, real_t(0));
      accr32.resize(n_cell, real_t(0));
      acnv20.resize(n_cell, real_t(0));
      acnv25.resize(n_cell, real_t(0));
      acnv32.resize(n_cell, real_t(0));
      // if(opts_init.adaptive_sstp_cond)
      //   sstp_cond_percell.resize(n_cell, opts_init.sstp_cond); // in adaptive substepping start with opts_init.sstp_cond substeps per cell
    }
  };
};
