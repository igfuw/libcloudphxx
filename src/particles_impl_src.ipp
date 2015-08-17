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
      template<typename t_a, typename t_b>
      struct two_keys_sort
      {
        BOOST_GPU_ENABLED
        bool operator()(const thrust::tuple<t_a, t_b> &a, const thrust::tuple<t_a, t_b> &b)
        {
          if(a.head < b.head) return true;
          if(a.head == b.head) return a.tail < b.tail;
          return false;
        }
      };
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

      // --- sort already existing SDs; primary key ijk, secondary rd ---
      // filling-in sorted_id with a sequence
      thrust::sequence(sorted_id.begin(), sorted_id.end());
     
      // use tmp_device_real_part as tmp copy of rd3
      thrust::copy(
        rd3.begin(), rd3.end(), // from
        tmp_device_real_part.begin()      // to
      );

      // copy ijk to sorted ijk
      thrust::copy(
        ijk.begin(), ijk.end(), // from
        sorted_ijk.begin()      // to
      );

      // sorting by ijk and rd3
      thrust::sort_by_key(
        thrust::make_zip_iterator(thrust::make_tuple(
          sorted_ijk.begin(),
          tmp_device_real_part.begin()
        )),
        thrust::make_zip_iterator(thrust::make_tuple(
          sorted_ijk.begin(),
          tmp_device_real_part.begin()
        )) + n_part,                           // keys
        sorted_id.begin(),                     // values
        detail::two_keys_sort<thrust_size_t, real_t>()
      ); 

      // analyze distribution to get rd_min and max needed for bin sizes
      // TODO: this could be done once at the beginning of the simulation
      dist_analysis(
        opts_init.src_dry_distros.begin()->second,
        opts_init.src_sd_conc,
        dt
      ); 
/*
      // --- see how many already existing SDs fall into size bins
      {
        // tmp vector with bin number of existing SDs
        thrust_device::vector<unsigned int> &bin_no(tmp_device_n_part);
 
        thrust::transform(
          rd3.begin(),
 
 
        // tmp vector to hold number of particles in a given size bin in a given cell
        // potentially could be rather large...
        // TODO: will it be emptied when it goes out of scope?
        thrust_device::vector<thrust_size_t> bin_cell_count(opts_init.src_sd_conc * n_cell);
        // TODO: some cells will have different bins due to differen vol, hense different src_sd_conc.....

        // tmp vector od IDs of SDs that are the smallest ones to fall into the first bin in given cell
        thrust_device::vector<thrust_size_t> &first_id(tmp_device_size_cell);
      }
*/
      // ------ update all vectors between n_part_old and n_part ------

      // init ijk and n_part, resize vectors
      // also set n_part_old and n_part_to_init used by init_dry and init_wet
      init_ijk();


      // init rd, n
      init_dry(
        opts_init.src_dry_distros.begin()->first,
        opts_init.src_dry_distros.begin()->second
      ); // TODO: document that n_of_lnrd_stp is expected!

      // init rw
      init_wet();

      // init x, y, z, i, j, k
      init_xyz();

      // init chem (TODO)
 
      // --- after source particles are no longer sorted ---
      sorted = false;
    }
  };  
};
