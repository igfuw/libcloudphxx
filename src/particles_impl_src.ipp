// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */
#include <limits>
#include <thrust/unique.h>

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

      template<typename real_t, typename n_t, typename res_t, int out_of_bins>
      struct get_bin_no
      {
        real_t log_rd_min, log_rd_max, log_diff;
        get_bin_no(const real_t &log_rd_min, const real_t &log_rd_max): log_rd_min(log_rd_min), log_rd_max(log_rd_max), log_diff(log_rd_max - log_rd_min) {}
        
        BOOST_GPU_ENABLED
        res_t operator()(real_t rd3, n_t count_num)
        {
          real_t log_rd = log(rd3) / 3.;
          if(log_rd < log_rd_min || log_rd >= log_rd_max || count_num == 0) return out_of_bins;
          else  return ( (log_rd - log_rd_min) / log_diff * count_num); //count num is the number of bins in the cell
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
      // TODO: do all of this only on SDs in cells below src_z1?
      // filling-in sorted_id with a sequence

      thrust::sequence(sorted_id.begin(), sorted_id.end());
     
      // tmp vector with sorted rd3
      thrust_device::vector<real_t> &sorted_rd3(tmp_device_real_part);

      // use sorted_rd3 as tmp copy of rd3
      thrust::copy(
        rd3.begin(), rd3.end(), // from
        sorted_rd3.begin()      // to
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
          sorted_rd3.begin()
        )),
        thrust::make_zip_iterator(thrust::make_tuple(
          sorted_ijk.begin(),
          sorted_rd3.begin()
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

      // --- see how many already existing SDs fall into size bins
      {
        namespace arg = thrust::placeholders;
        // tmp vector with bin number of existing SDs
        thrust_device::vector<unsigned int> &bin_no(tmp_device_n_part);

        enum {out_of_bins = std::numeric_limits<unsigned int>::max()};
        // calc bin no
        thrust::transform(
          sorted_rd3.begin(),
          sorted_rd3.end(),
          thrust::make_permutation_iterator(count_num.begin(), sorted_ijk.begin()),
          bin_no.begin(),
          detail::get_bin_no<real_t, n_t, unsigned int, out_of_bins>(log_rd_min, log_rd_max)
        );

        // set no of particles to init
        n_part_old = n_part;
        n_part_to_init = thrust::reduce(count_num.begin(), count_num.end());// - count_bins;

        // -- decrease count_num by number of SDs with increased n -- 
        {
//          debug::print(count_num);
          thrust_device::vector<unsigned int> tmp_bin_no(n_part);
          thrust::copy(bin_no.begin(), bin_no.end(), tmp_bin_no.begin());

          thrust_size_t n_out_of_bins = thrust::count(tmp_bin_no.begin(), tmp_bin_no.end(), out_of_bins);
        
//          printf("przed\n");
  //        debug::print(tmp_bin_no);
    //      debug::print(sorted_ijk);
          thrust::remove_if(
            thrust::make_zip_iterator(thrust::make_tuple(
              tmp_bin_no.begin(), 
              sorted_ijk.begin()
            )),
            thrust::make_zip_iterator(thrust::make_tuple(
              tmp_bin_no.begin(), 
              sorted_ijk.begin()
            )) + n_part,
            tmp_bin_no.begin(),
            arg::_1 == out_of_bins
          );
      //    printf("po usunieciu out_of_bins, n_out_of_bins = %lf\n", real_t(n_out_of_bins));
//          debug::print(tmp_bin_no);
//          debug::print(sorted_ijk);
          // sorted_ijk no longer valid

          thrust_size_t count_bins;
          {
            thrust::pair<
              thrust_device::vector<unsigned int>::iterator,
              typename thrust_device::vector<thrust_size_t>::iterator
            > n = thrust::unique_by_key(tmp_bin_no.begin(), tmp_bin_no.begin() + n_part - n_out_of_bins, sorted_ijk.begin());
            count_bins = n.first - tmp_bin_no.begin();
          }
//          printf("po unique, count_bins = %lf\n", real_t(count_bins));
//          debug::print(tmp_bin_no);
//          debug::print(sorted_ijk);

          thrust::fill(tmp_bin_no.begin(), tmp_bin_no.begin() + count_bins, 1);
          {
            thrust::pair<
              thrust_device::vector<thrust_size_t>::iterator,
              typename thrust_device::vector<unsigned int>::iterator
            > n =  thrust::reduce_by_key(
              sorted_ijk.begin(), sorted_ijk.begin() + count_bins, 
              tmp_bin_no.begin(), 
              sorted_ijk.begin(), 
              tmp_bin_no.begin()
            );
            count_bins = n.first - sorted_ijk.begin(); 
          }

//          printf("po fill i reduce, count_bins = %lf\n", real_t(count_bins));
//          debug::print(tmp_bin_no);
//          debug::print(sorted_ijk);

          thrust::transform(
            thrust::make_permutation_iterator(count_num.begin(), sorted_ijk.begin()),
            thrust::make_permutation_iterator(count_num.begin(), sorted_ijk.begin()) + count_bins,
            tmp_bin_no.begin(),
            thrust::make_permutation_iterator(count_num.begin(), sorted_ijk.begin()),
            arg::_1 - arg::_2
          );
//          debug::print(count_num);
        }

        // tmp vector to hold number of particles in a given size bin in a given cell
        // TODO: will it be emptied when it goes out of scope?
        thrust_device::vector<thrust_size_t> bin_cell_count(n_part_to_init +  n_cell + 1); // needs space for out_of_bins
        // tmp vector for number of particles in bins up to this one
        thrust_device::vector<thrust_size_t> bin_cell_count_ptr(n_part_to_init +  n_cell + 1);

        // calc no of SDs in bins/cells
        thrust::pair<
          thrust_device::vector<unsigned int>::iterator,
          typename thrust_device::vector<thrust_size_t>::iterator
        > n = thrust::reduce_by_key(
            bin_no.begin(),
            bin_no.begin() + n_part_old,
            thrust::make_constant_iterator<unsigned int>(1),
            bin_no.begin(),// output bin no - in place 
            bin_cell_count.begin()// output number of SDs
          );

        thrust_size_t count_bins = n.first - bin_no.begin(); // number of bins with SDs inside, includes the out_of_bins
        printf("po reduce - ile starych przypada na bin, count_bins = %lf\n", real_t(count_bins));
        debug::print(bin_no);
        debug::print(bin_cell_count);

        // number of SDs in bins up to (i-1)
        thrust::exclusive_scan(
          bin_cell_count.begin(), 
          bin_cell_count.begin() + count_bins, 
          bin_cell_count_ptr.begin()
        );

        // remove out of bins from bin_cell_count, bins_no
        count_bins = 
          thrust::remove_if(
             thrust::make_zip_iterator(thrust::make_tuple(
               bin_no.begin(), 
               bin_cell_count.begin()
             )),
             thrust::make_zip_iterator(thrust::make_tuple(
               bin_no.begin(), 
               bin_cell_count.begin()
             )) + count_bins,
             bin_no.begin(),
             arg::_1 == out_of_bins
          ) - 
           thrust::make_zip_iterator(thrust::make_tuple(
             bin_no.begin(), 
             bin_cell_count.begin()
           )); // count_bins now does not count out_of_bins

        printf("po usunieciu out_of_bins, count_bins = %lf\n", real_t(count_bins));
        debug::print(bin_no);
        debug::print(bin_cell_count);

        n_part_to_init -= count_bins;
        assert(n_part_to_init == thrust::reduce(count_num.begin(), count_num.end()));
        n_part += n_part_to_init;
        hskpng_resize_npart();   
        // dotad sie zgadza - nie sprawdzalem tylko ptr

        // randomly select which old SD will be increased
        // overwrites sorted_rd3
        rand_u01(count_bins);

        // TODO: merge the following transforms into one
        thrust::transform(
          u01.begin(),
          u01.begin() + count_bins,
          bin_cell_count.begin(),
          bin_cell_count.begin(),
          arg::_1 * arg::_2
        );
        thrust::transform(
          bin_cell_count_ptr.begin(),
          bin_cell_count_ptr.begin() + count_bins,
          bin_cell_count.begin(),
          bin_cell_count.begin(),
          arg::_1 + arg::_2
        );
        debug::print(bin_cell_count);

        // increase multiplicity of existing SDs
/*
        thrust::transform(
          thrust::make_permutation_iterator(
            n.begin(), thrust::make_permutation_iterator(
              sorted_id.begin(), bin_cell_count.begin()      // translates no of sorted SD into id
            )
          ),
          thrust::make_permutation_iterator(
            n.begin(), thrust::make_permutation_iterator(
              sorted_id.begin(), bin_cell_count.begin()
            )
          ) + count_bins,
          n.begin(),     // output
                         // operation
        );
*/
      // --- set rd3 of those that did not have counterparts ---



 
      }

      // ------ update all vectors between n_part_old and n_part ------

      // init ijk and n_part, resize vectors
      // also set n_part_old and n_part_to_init used by init_dry and init_wet
      init_ijk();

      // init n
      init_n(
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
