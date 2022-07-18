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
      eta.resize(n_cell.get()); 

      count_ijk.resize(n_cell.get());
      count_num.resize(n_cell.get());
      count_mom.resize(n_cell.get());
      count_n = 0;

      // initialising device temporary arrays
      tmp_device_real_cell.resize(n_cell.get());
      tmp_device_real_cell1.resize(n_cell.get());
      tmp_device_real_cell2.resize(n_cell.get());
      tmp_device_size_cell.resize(n_cell.get());
      tmp_host_size_cell.resize(n_cell.get());
      tmp_host_real_cell.resize(n_cell.get());
      if(allow_sstp_cond && !opts_init.exact_sstp_cond)
      {
        sstp_tmp_rv.resize(n_cell.get());
        sstp_tmp_th.resize(n_cell.get());
        sstp_tmp_rh.resize(n_cell.get());
      }

      // arrays on the refined grid
      if(opts_init.n_ref > 1)
      {
        T_ref.resize(n_cell_ref);
        RH_ref.resize(n_cell_ref); 
        eta_ref.resize(n_cell_ref); 
       
        tmp_device_real_cell_ref.resize(n_cell_ref);
        tmp_device_real_cell_ref1.resize(n_cell_ref);

        count_ijk_ref.resize(n_cell_ref);
        count_num_ref.resize(n_cell_ref);
        count_mom_ref.resize(n_cell_ref);
        count_n_ref = 0;
      }
    }
  };
};
