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
    void particles_t<real_t, device>::impl::fill_outbuf()
    {
      thrust::fill(tmp_host_real_cell.begin(), tmp_host_real_cell.end(), 0);

#if defined(__NVCC__)
      thrust::copy(
        count_ijk.begin(), count_ijk.end(), // from
        tmp_host_size_cell.begin()
      );
#endif

#if !defined(__NVCC__)
      thrust_device::vector<thrust_size_t> &pi(count_ijk);
#else
      thrust::host_vector<thrust_size_t> &pi(tmp_host_size_cell);
#endif

      thrust::copy(
	count_mom.begin(),               // input - begin
	count_mom.begin() + count_n,     // input - end
	thrust::make_permutation_iterator(  // output
	  tmp_host_real_cell.begin(),         // data
	  pi.begin()                          // permutation
	)
      );
    }

    template <typename real_t, backend_t device>
    std::vector<real_t> particles_t<real_t, device>::impl::fill_attr_outbuf(const std::string &name)
    {
      const std::set<std::string> attr_names = {"rw2", "rd3", "kappa", "rd3_insol", "T_freeze", "a_ice", "c_ice", "rho_i", "x", "y", "z"}; // TODO implement "n" - it is n_t type and others are real_t
      if (std::find(std::begin(attr_names), std::end(attr_names), name) == std::end(attr_names))
        throw std::runtime_error("Unknown attribute name passed to get_attr.");

      const thrust_device::vector<real_t> &dv(
        name == "rw2" ? rw2 : 
        name == "rd3" ? rd3 : 
        name == "kappa" ? kpa :
        name == "rd3_insol" ? rd3_insol :
        name == "T_freeze" ? T_freeze :
        name == "a_ice" ? a_ice :
        name == "c_ice" ? c_ice :
        name == "rho_i" ? rho_i :
        name == "x" ? x :
        name == "y" ? y :
        z); 

      // NOTE: for host backends (i.e. undefined __NVCC__) we could return the vector directly, without a copy;
      //       however, if output was done concurrently, values in the diagnosed vector might change after the call to fill_attr_outbuf.
      std::vector<real_t> out(n_part);
      thrust::copy(
        dv.begin(), dv.end(),
        out.begin()
      );
      return out;
    }
  };
};
