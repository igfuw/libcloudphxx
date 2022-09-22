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
 //     eta.resize(n_cell); 
 //     T.resize(n_cell); 
 //     RH.resize(n_cell); 

//      count_ijk.resize(n_cell);
//      count_num.resize(n_cell);
//      count_mom.resize(n_cell);
//      count_n.get() = 0;
//      count_n.get_ref() = 0;

      // initialising device temporary arrays
//      tmp_device_real_cell.resize(n_cell.get());
//      tmp_device_real_cell1.resize(n_cell);
//      tmp_device_real_cell2.resize(n_cell);
      tmp_device_size_cell.resize(n_cell.get());
//      tmp_host_size_cell.resize(n_cell.get());
//      tmp_host_real_cell.resize(n_cell.get());
      if(allow_sstp_cond && !opts_init.exact_sstp_cond)
      {
        sstp_tmp_rv.resize(0, n_cell.get_ref());
        sstp_tmp_th.resize(0, n_cell.get_ref());
        sstp_tmp_rh.resize(0, n_cell.get_ref());
      }

      //TODO: init ijk_ref2ijk. it's size is already initialized
    }
  };
};
