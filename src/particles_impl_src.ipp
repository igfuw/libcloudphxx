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
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::src(const real_t &dt)
    {   
      // sanity checks
      if(n_dims < 2)            throw std::runtime_error("Source only works in 2D and 3D");
      if(opts_init.chem_switch) throw std::runtime_error("Source is not yet compatible with chemistry.");

      // set number of SDs to init; use count_num as storage
      {
        namespace arg = thrust::placeholders;
        thrust::fill(count_num.begin(), count_num.end(), 0);
      
        thrust_size_t k1 = opts_init.src_z1 / opts_init.dz + 0.5; // k index of the heighest cell we create SDs in
        // some cells may be used only partially in thr super-droplet method
        // e.g. when Lagrangian domain (x0, x1, etc...) is smaller than the 
        // Eulerian domain (0, nx*dx, etc...)
        // sd_conc defines number of SDs per Eulerian cell
        thrust::transform_if(
          dv.begin(), dv.end(), 
          thrust::make_counting_iterator(0),
          count_num.begin(), 
          real_t(opts_init.src_sd_conc) *                           // no of SDs to create
          arg::_1 / (opts_init.dx * opts_init.dy * opts_init.dz) +  // ratio of volumes
          real_t(0.5),             
          (arg::_1 % opts_init.nz) < k1
        ); 
      }

      // TODO: assert that we do not introduce particles into supersaturated cells?

      // update all vectors between n_part_old and n_part
      // init x, y, z, ijk, i, j, k, n_part
      // also set n_part_old and n_part_to_init used by init_dry and init_wet
      init_xyz_helper();

      // init rd, n

      init_dry(
        opts_init.src_dry_distros.begin()->first,
        opts_init.src_dry_distros.begin()->second,
        opts_init.src_sd_conc,
        dt
      ); // TODO: document that n_of_lnrd_stp is expected!

      // init rw
      init_wet();

      // init chem (TODO)
 
      // after source particles are no longer sorted
      sorted = false;
    }
  };  
};
