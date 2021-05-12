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
    // create new aerosol particles to relax towards a size distribution
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::rlx_dry_distros()
    {   
      namespace arg = thrust::placeholders;

      // vectors of size nz used in calculation of horizontal averages, TODO: allocate them at init
      thrust::device_vector<real_t> hor_avg(opts_init.nz);
      thrust::device_vector<real_t> hor_avg_count(opts_init.nz);
      thrust::device_vector<int> hor_avg_k(opts_init.nz);
      thrust::device_vector<real_t> expected_hor_avg(opts_init.nz);

      // calc sum of ln(rd) ranges of all relax distributions
      real_t tot_lnrd_rng = 0.;
      for (typename opts_init_t<real_t>::rlx_dry_distros_t::const_iterator ddi = opts_init.rlx_dry_distros.begin(); ddi != opts_init.rlx_dry_distros.end(); ++ddi)
      {
        dist_analysis_sd_conc(
          *(std::get<0>(ddi->second)),
          opts_init.rlx_bins
        );
        tot_lnrd_rng += log_rd_max - log_rd_min;
      }

      // initialize SDs of each kappa-type
      for (typename opts_init_t<real_t>::rlx_dry_distros_t::const_iterator ddi = opts_init.rlx_dry_distros.begin(); ddi != opts_init.rlx_dry_distros.end(); ++ddi)
      {
        const auto &n_of_lnrd_stp(*(std::get<0>(ddi->second)));

        // analyze distribution to get rd_min and max needed for bin sizes
        // TODO: this was done a moment ago!
        // TODO2: this could be done once at the start of the simulation!
        dist_analysis_sd_conc(
          n_of_lnrd_stp,
          opts_init.rlx_bins
        );
        const real_t lnrd_rng = log_rd_max - log_rd_min;



        // calculate bin edges (in rd3)
        // tmp vector with bin edges, probably could be allocated once in init
        const int n_bins = opts_init.rlx_bins * lnrd_rng / tot_lnrd_rng;
        const real_t lnrd_bin_size = lnrd_rng / n_bins;
        std::vector<real_t> bin_rd3_left_edges(n_bins+1); // on CPU because of small number of edges
        std::iota(bin_rd3_left_edges.begin(), bin_rd3_left_edges.end(), 0); // fill with a 0,1,2,... sequence
        std::transform(bin_rd3_left_edges.begin(), bin_rd3_left_edges.end(), bin_rd3_left_edges.begin(), [log_rd_min_val=log_rd_min, lnrd_bin_size] (real_t bin_number) { return std::exp( 3 * (log_rd_min_val + bin_number * lnrd_bin_size)) ; }); // calculate left edges

        // minimum and maximum cell indices
        const int z_min_index = (std::get<2>(ddi->second)).first  / opts_init.nz + 0.5,
                  z_max_index = (std::get<2>(ddi->second)).second / opts_init.nz + 0.5;

        // loop over the bins
        for(int bin_number=0; bin_number<bin_rd3_left_edges.size()-1; ++bin_number)
        {
          const real_t rd3_min = bin_rd3_left_edges.at(bin_number),
                       rd3_max = bin_rd3_left_edges.at(bin_number+1);
          // TODO: these selections could be optimised
          // select droplets within the desired kappa range
          moms_rng((std::get<1>(ddi->second)).first, (std::get<1>(ddi->second)).second, kpa.begin(), false);
          // out of those, select droplets within the desired rd3 range
          moms_rng(rd3_min, rd3_max, rd3.begin(), true);
          // calculate 0-th non-specific moment of rd3 (number of droplets in a cell) of droplets in this rd3 and kappa range
          moms_calc(rd3.begin(), 0, false);
          // divide by volume
          thrust::transform(
            count_mom.begin(), count_mom.begin() + count_n,     // input - first arg
            thrust::make_permutation_iterator(                  // input - second arg
              dv.begin(),
              count_ijk.begin()
            ),
            count_mom.begin(),                                  // output (in place)
            thrust::divides<real_t>()
          );

          // horizontal average of this moment
          thrust::fill(hor_avg.begin(), hor_avg.end(), 0);
          thrust_device::vector<thrust_size_t> &count_k(tmp_device_size_cell);  // NOTE: tmp_device_size_cell is also used in some other inits, careful not to overwrite it!
          thrust::transform(count_ijk.begin(), count_ijk.begin() + count_n, count_k.begin(), arg::_1 % opts_init.nz);
          thrust::sort_by_key(count_k.begin(), count_k.begin() + count_n, count_mom.begin());
          auto new_end = thrust::reduce_by_key(count_k.begin(), count_k.begin() + count_n, count_mom.begin(), hor_avg_k.begin(), hor_avg_count.begin()); 
          int number_of_levels_with_droplets = new_end.first - hor_avg_k.begin();
          thrust::copy(hor_avg_count.begin(), hor_avg_count.begin() + number_of_levels_with_droplets, thrust::make_permutation_iterator(hor_avg.begin(), hor_avg_k.begin()));
          // divide sum by the number of cells at this level
          thrust::transform(hor_avg.begin(), hor_avg.end(), hor_avg.begin(), arg::_1 / (opts_init.nx * m1(opts_init.ny)));

          
          // calculate expected CCN conc
          const real_t bin_lnrd_center = log_rd_min + (bin_number + 0.5) * lnrd_bin_size;
          const real_t expected_STP_concentration = n_of_lnrd_stp(bin_lnrd_center) * lnrd_bin_size;
          thrust::fill(expected_hor_avg.begin(), expected_hor_avg.end(), expected_STP_concentration);
 
          // correcting STP -> actual ambient conditions
          if(!opts_init.aerosol_independent_of_rhod)
          {
            using common::earth::rho_stp;
            thrust::transform(
              expected_hor_avg.begin(), 
              expected_hor_avg.begin() + opts_init.nz, 
              rhod.begin(),                 // rhod has size ncell, but vertical cooridnate varies first, so rhod.begin() to rhod.begin()+nz should be the vertical profile?
              expected_hor_avg.begin(), 
              arg::_1 * arg::_2 / real_t(rho_stp<real_t>() / si::kilograms * si::cubic_metres)
            );
          }

          // set to zero outside of the defined range of altitudes
          thrust::transform_if(thrust::make_counting_iterator<int>(0), thrust::make_counting_iterator<int>(opts_init.nz), expected_hor_avg.begin(), arg::_1 = 0, arg::_1 < z_min_index || arg::_1 >= z_max_index); 

          std::cerr << "bin number: " << bin_number 
            << " rd_range: (" << std::pow(rd3_min, 1./3.) << ", " << std::pow(rd3_max, 1./3.) 
            << " r_center: " << std::exp(bin_lnrd_center) 
            << " z_indices: (" << z_min_index << ", " << z_max_index << "), " 
            << " expected STD concentration: " << expected_STP_concentration 
            << std::endl;
        
          std::cerr << "hor_avg:" << std::endl;
          debug::print(hor_avg);
        
          std::cerr << "expected_hor_avg:" << std::endl;
          debug::print(expected_hor_avg);

          // calculate how many CCN are missing
          thrust::device_vector<real_t> &hor_missing(expected_hor_avg);
          thrust::transform(expected_hor_avg.begin(), expected_hor_avg.end(), hor_avg.begin(), hor_missing.begin(), arg::_1 - arg::_2);
          thrust::replace_if(hor_missing.begin(), hor_missing.end(), arg::_1 < 0, 0);
         
          std::cerr << "hor_missing:" << std::endl;
          debug::print(hor_missing);

          // TODO: watch out not to mess up sorting while adding SDs to the bins, because moms_X functions require sorted data...

        }
      }
        

        // 

//      // set number of SDs to init; use count_num as storage
//      init_count_num_src(opts_init.src_sd_conc);

/*

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
            > it_pair =  thrust::reduce_by_key(
              sorted_ijk.begin(), sorted_ijk.begin() + count_bins, 
              thrust::make_constant_iterator<thrust_size_t>(1), 
              sorted_ijk.begin(), 
              tmp_bin_no.begin()
            );
            count_bins = it_pair.first - sorted_ijk.begin(); // now it counts no of cells that have any bins matched
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
      */
    }
  };  
};
