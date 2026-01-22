// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */
// #include <limits>
#include <thrust/unique.h>
#include <thrust/binary_search.h>
#include <numeric>


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

      template<class real_t>
      struct calc_n_sd_to_create
      {
        const real_t tolerance;
        const int n_sd_per_bin;

        calc_n_sd_to_create(const real_t tol, const real_t n_sd_per_bin_real):
          tolerance(tol),
          n_sd_per_bin(std::max<int>(1, int(n_sd_per_bin_real + 0.5)))
          {}

        BOOST_GPU_ENABLED
        int operator()(const real_t &a1, const real_t &a2)
        {
          return a2 > real_t(0) ?
            a1 / a2 > tolerance ? n_sd_per_bin : 0
          : 0;
        }
      };

      // domain volume at this height level
      template <typename real_t>
      struct hor_dv_eval
      {
        // note: having a copy of opts_init here causes CUDA crashes (alignment problems?)
        const real_t
          dz,
          x0, y0, z0,
          x1, y1, z1;

        hor_dv_eval(const opts_init_t<real_t> &o) :
          dz(o.dz),
          x0(o.x0), y0(o.y0), z0(o.z0),
          x1(o.x1), y1(o.y1), z1(o.z1)
        {}

        BOOST_GPU_ENABLED
        real_t operator()(const int &k)
        {
#if !defined(__NVCC__)
          using std::min;
          using std::max;
#endif
          return
            max(real_t(0),
              (x1 - x0) *
              (y1 - y0) * // NOTE: size in y is taken into account even in 2D!
              (min((k + 1) * dz, z1) - max(k * dz, z0))
            );
        }
      };
    };

    // create new aerosol particles to relax towards a size distribution
    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::rlx_dry_distros(const real_t dt)
    {   
      namespace arg = thrust::placeholders;

      // vectors of size nz used in calculation of horizontal averages, TODO: allocate them at init
      thrust_device::vector<real_t> hor_sum(opts_init.nz);
      thrust_device::vector<real_t> hor_sum_count(opts_init.nz);
      thrust_device::vector<real_t> hor_missing(opts_init.nz);
      thrust_device::vector<thrust_size_t> hor_sum_k(opts_init.nz);
      thrust_device::vector<real_t> expected_hor_sum(opts_init.nz);
      thrust_device::vector<thrust_size_t> n_SD_to_create(opts_init.nz); // could be bool, but then thrust::reduce does not add bools as expected

      // calc sum of ln(rd) ranges of all relax distributions
      real_t tot_lnrd_rng = 0.;
      for (typename opts_init_t<real_t>::rlx_dry_distros_t::const_iterator ddi = opts_init.rlx_dry_distros.begin(); ddi != opts_init.rlx_dry_distros.end(); ++ddi)
      {
        init_dist_analysis_sd_conc(
          *(std::get<0>(ddi->second)),
          opts_init.rlx_bins
        );
        tot_lnrd_rng += log_rd_max - log_rd_min;
      }

      const auto n_part_pre_relax = n_part;

      // initialize SDs of each kappa-type
      for (typename opts_init_t<real_t>::rlx_dry_distros_t::const_iterator ddi = opts_init.rlx_dry_distros.begin(); ddi != opts_init.rlx_dry_distros.end(); ++ddi)
      {
        const auto &kappa(ddi->first);
        assert(kappa >= 0);
        const auto &n_of_lnrd_stp(*(std::get<0>(ddi->second)));

        // analyze distribution to get rd_min and max needed for bin sizes
        // TODO: this was done a moment ago!
        // TODO2: this could be done once at the start of the simulation!
        init_dist_analysis_sd_conc(
          n_of_lnrd_stp,
          opts_init.rlx_bins
        );
        const real_t lnrd_rng = log_rd_max - log_rd_min;
        assert(lnrd_rng > 0);

        // calculate bin edges (in rd3)
        // tmp vector with bin edges, probably could be allocated once in init
        const int n_bins = opts_init.rlx_bins * lnrd_rng / tot_lnrd_rng;
        assert(n_bins>0);
        const real_t lnrd_bin_size = lnrd_rng / n_bins;
        assert(lnrd_bin_size > 0);
        std::vector<real_t> bin_rd3_left_edges(n_bins+1); // on CPU because of small number of edges
        std::iota(bin_rd3_left_edges.begin(), bin_rd3_left_edges.end(), 0); // fill with a 0,1,2,... sequence
        std::transform(bin_rd3_left_edges.begin(), bin_rd3_left_edges.end(), bin_rd3_left_edges.begin(), [log_rd_min_val=log_rd_min, lnrd_bin_size] (real_t bin_number) { return std::exp( 3 * (log_rd_min_val + bin_number * lnrd_bin_size)) ; }); // calculate left edges

        // minimum and maximum cell indices
        const int z_min_index = (std::get<2>(ddi->second)).first  / opts_init.dz,
                  z_max_index = (std::get<2>(ddi->second)).second / opts_init.dz;

        assert(z_max_index >= z_min_index);
        assert(z_min_index >= 0);
        assert(z_max_index < opts_init.nz);

        const auto n_part_pre_bins_loop = n_part;

        real_t expected_STP_concentration_tot = 0;

        // loop over the bins
        for(int bin_number=0; bin_number<bin_rd3_left_edges.size()-1; ++bin_number)
        {
          const real_t rd3_min = bin_rd3_left_edges.at(bin_number),
                       rd3_max = bin_rd3_left_edges.at(bin_number+1);
          assert(rd3_min < rd3_max);
          // TODO: these selections could be optimised
          // select droplets within the desired kappa range; only from the droplets existing before relaxation began, because relaxed ones are not sorted
          moms_rng((std::get<1>(ddi->second)).first, (std::get<1>(ddi->second)).second, kpa.begin(), n_part_pre_relax, false);
          // out of those, select droplets within the desired rd3 range
          moms_rng(rd3_min, rd3_max, rd3.begin(), n_part_pre_relax, true);
          // calculate 0-th non-specific moment of rd3 (number of droplets in a cell) of droplets in this rd3 and kappa range
          moms_calc(rd3.begin(), n_part_pre_relax, 0, false);
          n_filtered_gp.reset(); // n_filtered not needed anymore

          // horizontal sum of this moment
          thrust::fill(hor_sum.begin(), hor_sum.end(), 0);
          {
            auto count_k_g = tmp_device_size_cell.get_guard();
            thrust_device::vector<thrust_size_t> &count_k(count_k_g.get());
            thrust::transform(count_ijk.begin(), count_ijk.begin() + count_n, count_k.begin(), arg::_1 % opts_init.nz);
            thrust::sort_by_key(count_k.begin(), count_k.begin() + count_n, count_mom.begin());

            auto new_end = thrust::reduce_by_key(count_k.begin(), count_k.begin() + count_n, count_mom.begin(), hor_sum_k.begin(), hor_sum_count.begin()); 

            int number_of_levels_with_droplets = new_end.first - hor_sum_k.begin(); // number of levels with any SD, not with SD in this size and kappa range
            
            assert(number_of_levels_with_droplets <= opts_init.nz);
            thrust::copy(hor_sum_count.begin(), hor_sum_count.begin() + number_of_levels_with_droplets, thrust::make_permutation_iterator(hor_sum.begin(), hor_sum_k.begin()));
          }
          // divide sum by the number of cells at this level
//          thrust::transform(hor_sum.begin(), hor_sum.end(), hor_sum.begin(), arg::_1 / (opts_init.nx * m1(opts_init.ny)));
          
          // calculate expected CCN number
          const real_t bin_lnrd_center = log_rd_min + (bin_number + 0.5) * lnrd_bin_size;
          const real_t expected_STP_concentration = n_of_lnrd_stp(bin_lnrd_center) * lnrd_bin_size;
          assert(expected_STP_concentration >= 0);
          expected_STP_concentration_tot += expected_STP_concentration;
          thrust::transform(zero, zero + opts_init.nz, expected_hor_sum.begin(), detail::hor_dv_eval<real_t>(opts_init)); // fill with volume of the domain at this level
          thrust::transform(expected_hor_sum.begin(), expected_hor_sum.end(), expected_hor_sum.begin(), expected_STP_concentration * arg::_1); // multiply by the expected concentration

          // TODO: check for overflows?
 
          // correcting STP -> actual ambient conditions
          if(!opts_init.aerosol_independent_of_rhod)
          {
            using common::earth::rho_stp;
            thrust::transform(
              expected_hor_sum.begin(), 
              expected_hor_sum.begin() + opts_init.nz, 
              rhod.begin(),                 // rhod has size ncell, but vertical cooridnate varies first, so rhod.begin() to rhod.begin()+nz should be the vertical profile?
              expected_hor_sum.begin(), 
              arg::_1 * arg::_2 / real_t(rho_stp<real_t>() / si::kilograms * si::cubic_metres)
            );
          }

          // set to zero outside of the defined range of altitudes
          thrust::replace_if(expected_hor_sum.begin(), expected_hor_sum.begin()+opts_init.nz, zero, arg::_1 < z_min_index || arg::_1 >= z_max_index, real_t(0));

          // calculate how many CCN are missing
          thrust::transform(expected_hor_sum.begin(), expected_hor_sum.end(), hor_sum.begin(), hor_missing.begin(), arg::_1 - arg::_2);
          thrust::replace_if(hor_missing.begin(), hor_missing.end(), arg::_1 < 0, 0);
         
          // set number of SDs to init; create only if concentration is lower than expected with a tolerance
          thrust::transform(hor_missing.begin(), hor_missing.end(), expected_hor_sum.begin(), n_SD_to_create.begin(), detail::calc_n_sd_to_create<real_t>(config.rlx_conc_tolerance, opts_init.rlx_sd_per_bin));

          n_part_old = n_part;
          n_part_to_init = thrust::reduce(n_SD_to_create.begin(), n_SD_to_create.end());
          n_part = n_part_old + n_part_to_init;

          // resize arrays set in the bins loop: cell indices and rd3, resize should be cheap, because we allocate a large chunk of memory at the start
          hskpng_resize_npart();

          // --- init k ---
          {
            auto ptr_g = tmp_device_size_cell.get_guard();
            thrust_device::vector<thrust_size_t> &ptr(ptr_g.get());
            thrust::exclusive_scan(n_SD_to_create.begin(), n_SD_to_create.end(), ptr.begin()); // number of SDs in cells to init up to (i-1)

            thrust_device::vector<thrust_size_t> &k(k_gp->get());
            thrust::for_each(
              thrust::make_zip_iterator(thrust::make_tuple(
                n_SD_to_create.begin(), ptr.begin(), zero
              )),
              thrust::make_zip_iterator(thrust::make_tuple(
                n_SD_to_create.end(), ptr.end(), zero + opts_init.nz
              )),
              detail::arbitrary_sequence<thrust_size_t>(&(k[n_part_old]))
            );


            // --- init multiplicities (includes casting from real to n) ---
#if !defined(__NVCC__)
            using std::min;
#endif

            thrust::for_each(
              thrust::make_zip_iterator(thrust::make_tuple(
                n_SD_to_create.begin(), ptr.begin(), 
                thrust::make_transform_iterator(hor_missing.begin(), arg::_1 / real_t(opts_init.rlx_sd_per_bin) * min(dt / opts_init.rlx_timescale, real_t(1)) + real_t(0.5))
              )),
              thrust::make_zip_iterator(thrust::make_tuple(
                n_SD_to_create.end(), ptr.end(),
                thrust::make_transform_iterator(hor_missing.end(), arg::_1 / real_t(opts_init.rlx_sd_per_bin) * min(dt / opts_init.rlx_timescale, real_t(1)) + real_t(0.5))
              )),
              detail::arbitrary_sequence<n_t>(&(n[n_part_old]))
            );
          }

          // detecting possible overflows of n type
          {
            thrust_size_t ix = thrust::max_element(n.begin() + n_part_old, n.end()) - (n.begin() + n_part_old);
            assert(n[ix] < (typename impl::n_t)(-1) / 10000);
          }

          // --- init of i and j ---
          // tossing random numbers [0,1)  TODO: do it once for all bins
          auto u01g = tmp_device_real_part.get_guard();
          thrust_device::vector<real_t> &u01 = u01g.get();
          rand_u01(u01, n_part_to_init * (n_dims)); // random numbers for: i, rd, j (j only in 3D)

          thrust_device::vector<thrust_size_t> &i(i_gp->get());
          thrust::transform(u01.begin(), u01.begin() + n_part_to_init, i.begin() + n_part_old, detail::multiply_by_constant_and_cast<real_t, thrust_size_t>(opts_init.nx));

          if(n_dims==3)
          {
            thrust_device::vector<thrust_size_t> &j(j_gp->get());
            thrust::transform(u01.begin() + 2*n_part_to_init, u01.begin() + 3*n_part_to_init, j.begin() + n_part_old, detail::multiply_by_constant_and_cast<real_t, thrust_size_t>(opts_init.ny));
          }

          // raveling i, j & k into ijk; only of the new SD
          ravel_ijk(n_part_old);

          // init dry radius
          // set rd3 randomized within the bin, uniformly distributed on the log(rd) axis
          const real_t lnrd_min = std::log(std::pow(rd3_min, 1./3.));
          thrust::transform(u01.begin() + 1*n_part_to_init, u01.begin() + 2*n_part_to_init, rd3.begin()+n_part_old, lnrd_min + arg::_1 * lnrd_bin_size);

          // converting from lnrd to rd3
          thrust::transform(
            rd3.begin() + n_part_old,
            rd3.end(),
            rd3.begin() + n_part_old,
            detail::exp3x<real_t>()
          );

          // NOTE: watch out not to mess up sorting while adding SDs to the bins, because moms_X functions require sorted data...
        } // end of the bins loop

        // init other SD characteristics that don't have to be initialized in the bins loop
        n_part_old = n_part_pre_bins_loop;
        n_part_to_init = n_part - n_part_old;
//        hskpng_resize_npart();

        init_SD_with_distros_finalize(lgrngn::kappa_rd_insol_t<real_t>{kappa, real_t(0)}, false); // no need to unravel ijk there, because i j k are already initialized
                                                     // TODO: we assume that relaxation produces water (not ice)

        // TODO: asserts of newly added SD parameters? e.g. how many SD, how big is multiplicity etc.
      } // end of the distros loop
      sorted = false;
    }
  };  
};
