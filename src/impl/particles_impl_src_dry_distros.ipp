// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */
// #include <limits>
#include <thrust/unique.h>
#include <thrust/binary_search.h>

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

      template<typename real_t, typename n_t, typename res_t, n_t out_of_bins>
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

    // create new aerosol particles based on a size distribution
    // if any SDs with dry radius similar to the one to be added are present,
    // we increase their multiplicity instead of adding new SDs
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::src_dry_distros(const real_t &dt)
    {   
      // set number of SDs to init; use count_num as storage
      init_count_num_src(opts_init.src_sd_conc);

      // -------- TODO: match not only sizes of old particles, but also kappas and chemical composition... --------

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
      dist_analysis_sd_conc(
        *(opts_init.src_dry_distros.begin()->second),
        opts_init.src_sd_conc,
        dt
      ); 

      // --- see how many already existing SDs match size bins ---
      {
        namespace arg = thrust::placeholders;

        // set no of particles to init
        n_part_old = n_part;
        n_part_to_init = thrust::reduce(count_num.begin(), count_num.end());
        n_part = n_part_old + n_part_to_init;
        hskpng_resize_npart();

        thrust_size_t n_part_bfr_src = n_part_old,
                      n_part_tot_in_src = n_part_to_init;

        // tmp vector with bin number of existing SDs
        thrust_device::vector<thrust_size_t> bin_no(n_part);

        const thrust_size_t out_of_bins = 4444444444; // would cause an error for src_sd_conc > out_of_bins
        // calc bin no
        thrust::transform(
          sorted_rd3.begin(),
          sorted_rd3.end(),
          thrust::make_permutation_iterator(count_num.begin(), sorted_ijk.begin()),
          bin_no.begin(),
          detail::get_bin_no<real_t, n_t, thrust_size_t, out_of_bins>(log_rd_min, log_rd_max)
        );
                    
        // init ijk and rd3 of new particles
        init_ijk();
        init_dry_sd_conc(); 

        // translate these new rd3 into bin no; bin_no just got resized
        thrust::transform(
          rd3.begin() + n_part_old,
          rd3.end(),
          thrust::make_permutation_iterator(count_num.begin(), ijk.begin() + n_part_old),
          bin_no.begin() + n_part_old,
          detail::get_bin_no<real_t, n_t, thrust_size_t, out_of_bins>(log_rd_min, log_rd_max)
        );

        // -- init new SDs that didnt have a match -- 
        {
          thrust_device::vector<thrust_size_t> tmp_bin_no(n_part_old);
          thrust::copy(bin_no.begin(), bin_no.begin() + n_part_old, tmp_bin_no.begin());

          thrust_size_t n_out_of_bins = thrust::count(tmp_bin_no.begin(), tmp_bin_no.end(), out_of_bins);
        
          // remove reference to those outside of bins from tmp_bin_no and sorted_ijk
          thrust::remove_if(
            sorted_ijk.begin(),
            sorted_ijk.begin() + n_part_old,
            tmp_bin_no.begin(),
            arg::_1 == out_of_bins
          );

          thrust::remove(
            tmp_bin_no.begin(),
            tmp_bin_no.begin() + n_part_old,
            out_of_bins
          ); // if these two removes are done in a single step with a tuple, it fails on CUDA; TODO: report this?

          thrust_size_t count_bins;
          {
          // remove duplicates from tmp_bin_no
            thrust::pair<
              thrust_device::vector<thrust_size_t>::iterator,
              typename thrust_device::vector<thrust_size_t>::iterator
            > np = thrust::unique_by_key(tmp_bin_no.begin(), tmp_bin_no.begin() + n_part_old - n_out_of_bins, sorted_ijk.begin());
            count_bins = np.first - tmp_bin_no.begin(); // total no of bins with a match
          }

          // --- remove rd3 and ijk of newly added SDs that have counterparts ---
          thrust_device::vector<bool> have_match(n_part_to_init);
          // find those with a match
          thrust::binary_search(
            thrust::make_zip_iterator(thrust::make_tuple(
              sorted_ijk.begin(),
              tmp_bin_no.begin()
            )),
            thrust::make_zip_iterator(thrust::make_tuple(
              sorted_ijk.begin(),
              tmp_bin_no.begin()
            )) + count_bins,
            thrust::make_zip_iterator(thrust::make_tuple(
              ijk.begin() + n_part_old,
              bin_no.begin() + n_part_old
            )),
            thrust::make_zip_iterator(thrust::make_tuple(
              ijk.begin() + n_part_old,
              bin_no.begin() + n_part_old
            )) + n_part_to_init,
            have_match.begin(),
            detail::two_keys_sort<thrust_size_t, thrust_size_t>()
          );
          // remove those with a match
          thrust::remove_if(
            thrust::make_zip_iterator(thrust::make_tuple(
              rd3.begin() + n_part_old,
              ijk.begin() + n_part_old
            )),
            thrust::make_zip_iterator(thrust::make_tuple(
              rd3.begin() + n_part_old,
              ijk.begin() + n_part_old
            )) + n_part_to_init,
            have_match.begin(),
            thrust::identity<bool>()
          );

          n_part_to_init -= count_bins;
          n_part -= count_bins;
          hskpng_resize_npart();

          // init other peoperties of SDs that didnt have a match
          init_kappa(
            opts_init.src_dry_distros.begin()->first
          ); 

          if(opts_init.diag_incloud_time)
            init_incloud_time();

          init_n_sd_conc(
            *(opts_init.src_dry_distros.begin()->second)
          ); // TODO: document that n_of_lnrd_stp is expected!

          // init rw
          init_wet();

          // init x, y, z, i, j, k
          init_xyz();

          // TODO: init chem
            
          {
            // count number of matched bins per cell
            thrust::pair<
              thrust_device::vector<thrust_size_t>::iterator,
              typename thrust_device::vector<thrust_size_t>::iterator
            > np =  thrust::reduce_by_key(
              sorted_ijk.begin(), sorted_ijk.begin() + count_bins, 
              thrust::make_constant_iterator<thrust_size_t>(1), 
              sorted_ijk.begin(), 
              tmp_bin_no.begin()
            );
            count_bins = np.first - sorted_ijk.begin(); // now it counts no of cells that have any bins matched
          }

          // set count_num to the number of SDs matched per cell
          // they still need to be initialized
          thrust::copy(
            tmp_bin_no.begin(),
            tmp_bin_no.begin() + count_bins,
            thrust::make_permutation_iterator(count_num.begin(), sorted_ijk.begin())
          );
          // sorted_ijk no longer valid
        }

        // tmp vector to hold number of particles in a given size bin in a given cell
        thrust_device::vector<thrust_size_t> bin_cell_count(n_part_tot_in_src +  n_cell + 1); // needs space for out_of_bins
        // tmp vector for number of particles in bins up to this one
        thrust_device::vector<thrust_size_t> bin_cell_count_ptr(n_part_tot_in_src +  n_cell + 1);

        thrust_size_t count_bins;
        {
          thrust_device::vector<thrust_size_t> &out(bin_cell_count_ptr); // use it temporarily
          // calc no of SDs in bins/cells
          thrust::pair<
            thrust_device::vector<thrust_size_t>::iterator,
            typename thrust_device::vector<thrust_size_t>::iterator
          > np = thrust::reduce_by_key(
              bin_no.begin(),
              bin_no.begin() + n_part_bfr_src,
              thrust::make_constant_iterator<thrust_size_t>(1),
              out.begin(),// output bin no - in place didn't work well, why?
              bin_cell_count.begin()// output number of SDs
            );
          count_bins = np.second - bin_cell_count.begin(); // number of bins with SDs inside, includes the out_of_bins
          thrust::copy(out.begin(), out.begin() + count_bins, bin_no.begin());
        }

        // number of SDs (incl. out_of_bins) in bins up to (i-1)
        thrust::exclusive_scan(
          bin_cell_count.begin(), 
          bin_cell_count.begin() + count_bins, 
          bin_cell_count_ptr.begin()
        );

        // remove out of bins from bin_cell_count, bins_no and bin cell_count_ptr
        count_bins = 
          thrust::remove_if(
             thrust::make_zip_iterator(thrust::make_tuple(
               bin_no.begin(), 
               bin_cell_count.begin(),
               bin_cell_count_ptr.begin()
             )),
             thrust::make_zip_iterator(thrust::make_tuple(
               bin_no.begin(), 
               bin_cell_count.begin(),
               bin_cell_count_ptr.begin()
             )) + count_bins,
             bin_no.begin(),
             arg::_1 == out_of_bins
          ) - 
           thrust::make_zip_iterator(thrust::make_tuple(
             bin_no.begin(), 
             bin_cell_count.begin(),
             bin_cell_count_ptr.begin()
           )); // count_bins now does not count out_of_bins

        // randomly select which old SD will be increased
        // overwrites sorted_rd3
        rand_u01(count_bins);

        // TODO: merge the following transforms into one
        
        // randomly choose one SD per bin
        thrust::transform(
          u01.begin(),
          u01.begin() + count_bins,
          bin_cell_count.begin(),
          bin_cell_count.begin(),
          thrust::multiplies<real_t>()
        );
        // translate no in bin to the total id
        thrust::transform(
          bin_cell_count_ptr.begin(),
          bin_cell_count_ptr.begin() + count_bins,
          bin_cell_count.begin(),
          bin_cell_count.begin(),
          thrust::plus<real_t>()
        );

        // --- increase multiplicity of existing SDs ---

        n_part_old = n_part; // number of those before src + no of those w/o match
        n_part_to_init = count_bins; // number of matched SDs
        n_part = n_part_old + n_part_to_init;
        hskpng_resize_npart();

        // copy rd3 and ijk of the selected SDs to the end of the respective vectors
        thrust::copy(
          thrust::make_permutation_iterator(
            rd3.begin(), thrust::make_permutation_iterator(
              sorted_id.begin(), bin_cell_count.begin()      // translates no of sorted SD into id
            )
          ),
          thrust::make_permutation_iterator(
            rd3.begin(), thrust::make_permutation_iterator(
              sorted_id.begin(), bin_cell_count.begin()
            )
          ) + count_bins,
          rd3.begin() + n_part_old     // output
        );
        thrust::copy(
          thrust::make_permutation_iterator(
            ijk.begin(), thrust::make_permutation_iterator(
              sorted_id.begin(), bin_cell_count.begin()      // translates no of sorted SD into id
            )
          ),
          thrust::make_permutation_iterator(
            ijk.begin(), thrust::make_permutation_iterator(
              sorted_id.begin(), bin_cell_count.begin()
            )
          ) + count_bins,
          ijk.begin() + n_part_old     // output
        );

        // init n of the copied SDs, but using the src distribution
        init_n_sd_conc(
          *(opts_init.src_dry_distros.begin()->second)
        ); // TODO: document that n_of_lnrd_stp is expected!

        // add the just-initialized multiplicities to the old ones
        thrust::transform(
          n.begin() + n_part_old,
          n.end(),
          thrust::make_permutation_iterator(
            n.begin(), thrust::make_permutation_iterator(
              sorted_id.begin(), bin_cell_count.begin()      // translates no of sorted SD into id
            )
          ),
          thrust::make_permutation_iterator(
            n.begin(), thrust::make_permutation_iterator(
              sorted_id.begin(), bin_cell_count.begin()      // translates no of sorted SD into id
            )
          ), //in-place
          thrust::plus<n_t>()
        );
        // TODO: check for overflows of na after addition

        // --- properly reduce size of the vectors back to no before src + no w/o match ---
        n_part = n_part_old;
        hskpng_resize_npart();   
      }
    }
  };  
};
