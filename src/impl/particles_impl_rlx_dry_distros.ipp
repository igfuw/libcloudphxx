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
      /// @brief returns ret_t(x*c)
      template <typename arg_t, typename ret_t>
      struct multiply_by_constant_and_cast
      {
        arg_t c;
        multiply_by_constant_and_cast(arg_t c) : c(c) {}

        BOOST_GPU_ENABLED
        ret_t operator()(arg_t x)
        {
          return ret_t(x*c);
        }
      };
    };

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
      thrust::device_vector<bool> create_SD(opts_init.nz);

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
        const int z_min_index = (std::get<2>(ddi->second)).first  / opts_init.nz,
                  z_max_index = (std::get<2>(ddi->second)).second / opts_init.nz;

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
          thrust::transform_if(zero, zero+opts_init.nz, expected_hor_avg.begin(), arg::_1 = 0, arg::_1 < z_min_index || arg::_1 > z_max_index); 

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
        
          // set number of SDs to init; create only if concentration is lower than expected with a tolerance
          thrust::transform(hor_missing.begin(), hor_missing.end(), expected_hor_avg.begin(), create_SD.begin(), arg::_2 > 0 && arg::_1 / arg::_2 > config.rlx_conc_tolerance); // WARNING: watch out for div by 0
         
          std::cerr << "create_SD:" << std::endl;
          debug::print(create_SD);

          n_part_old = n_part;
          n_part_to_init = thrust::reduce(create_SD.begin(), create_SD.end());
          n_part = n_part_old + n_part_to_init;

          // resize cell indices, resize should be cheap, because we allocate a large chunk of memory at the start
          ijk.resize(n_part);
          i.resize(n_part);
          k.resize(n_part);
          if(n_dims==3) j.resize(n_part); // we dont check in i and k because relax works only in 2D and 3D

          // k index based on create_SD
          thrust::copy_if(zero, zero+opts_init.nz, create_SD.begin(), k.begin()+n_part_old, arg::_1 == 1);
          // i and j random 
          // tossing random numbers [0,1)  TODO: do it once for all bins
          rand_u01(n_part_to_init * (n_dims-1));

          std::cerr << "u01:" << std::endl;
          debug::print(u01.begin(), u01.begin()+n_part_to_init);

          thrust::transform(u01.begin(), u01.begin() + n_part_to_init, i.begin() + n_part_old, detail::multiply_by_constant_and_cast<real_t, thrust_size_t>(opts_init.nx));
          if(n_dims==3) thrust::transform(u01.begin() + n_part_to_init, u01.begin() + 2*n_part_to_init, j.begin() + n_part_old, detail::multiply_by_constant_and_cast<real_t, thrust_size_t>(opts_init.ny));

          // raveling i, j & k into ijk; only of the new SD
          ravel_ijk(n_part_old);

          // set count_num to the number of SD to init per cell
          thrust::fill(count_num.begin(), count_num.end(), 0);
          thrust::scatter(thrust::make_constant_iterator<n_t>(1), thrust::make_constant_iterator<n_t>(1) + n_part_to_init, ijk.begin() + n_part_old, count_num.begin());

          std::cerr << "i:" << std::endl;
          debug::print(i.begin()+n_part_old, i.end());

          std::cerr << "k:" << std::endl;
          debug::print(k.begin()+n_part_old, k.end());

          std::cerr << "ijk:" << std::endl;
          debug::print(ijk.begin()+n_part_old, ijk.end());

          std::cerr << "count_num:" << std::endl;
          debug::print(count_num);

          // init i,j,k, (random i and j, k from create SD)
          // init count num
          // init ijk
          // other inits, TODO: how to make these init not in the bins loop?

          // TODO: watch out not to mess up sorting while adding SDs to the bins, because moms_X functions require sorted data...

          // at the end we need to set sorting=false

        }
      }

    }
  };  
};
