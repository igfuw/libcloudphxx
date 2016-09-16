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
      enum{na_ge_nb = -2, nb_gt_na = -1};

      struct summator
      {
        template<class tpl_t>
        BOOST_GPU_ENABLED
        void operator()(tpl_t tpl)
        {
          if(thrust::get<2>(tpl) <= 0) return; // do nothing if no collisions or first one passed was a SD with an uneven number in the cell
          thrust::get<3>(tpl) == na_ge_nb ?    // does the first SD of the pair have greater multiplicity?
            thrust::get<1>(tpl) += thrust::get<2>(tpl) * thrust::get<0>(tpl): // add col_no *val(SD with greater multiplicity) to the one with smaller multiplicity
            thrust::get<0>(tpl) += thrust::get<2>(tpl) * thrust::get<1>(tpl); // add col_no *val(SD with greater multiplicity) to the one with smaller multiplicity
        }
      };

      template<class real_t>
      struct weighted_summator
      {
        template<class tpl_t>
        BOOST_GPU_ENABLED
        void operator()(tpl_t tpl)
        {
          if(thrust::get<4>(tpl) <= 0) return; // do nothing if no collisions or first one passed was a SD with an uneven number in the cell

          // previous value of dry radius of the one with smaller multiplicity, it was allready updated in collide
          real_t rd3_old = 
            thrust::get<5>(tpl) == na_ge_nb ?    // does the first SD of the pair have greater multiplicity?
              thrust::get<3>(tpl) - thrust::get<4>(tpl) * thrust::get<2>(tpl) : // rd3_old = rd3_new - col_no * rd3_old_a
              thrust::get<2>(tpl) - thrust::get<4>(tpl) * thrust::get<3>(tpl);  // rd3_old = rd3_new - col_no * rd3_old_b

          for(int ci=0; ci<thrust::get<4>(tpl); ++ci) // loop over collisions between the pair
          {
            if(thrust::get<5>(tpl) == na_ge_nb)    // first SD of the pair has greater multiplicity
            {
              // update kappa_b
              thrust::get<1>(tpl) =  (thrust::get<0>(tpl) * thrust::get<2>(tpl) + // kappa_a * rd3_a
                                     thrust::get<1>(tpl) * rd3_old) /             // kappa_b * rd3_b_previous
                                     (thrust::get<2>(tpl) + rd3_old);             // rd3_a + rd3_b_previous
              // update rd3_old after collision
              rd3_old += thrust::get<2>(tpl); // add rd3_a
            }
            else    // second SD of the pair has greater multiplicity
            {
              // update kappa_a
              thrust::get<0>(tpl) =  (thrust::get<1>(tpl) * thrust::get<3>(tpl) + // kappa_b * rd3_b
                                     thrust::get<0>(tpl) * rd3_old) /             // kappa_a * rd3_a_previous
                                     (thrust::get<3>(tpl) + rd3_old);             // rd3_b + rd3_b_previous
              // update rd3_old after collision
              rd3_old += thrust::get<3>(tpl); // add rd3_b
            }
          }
        }
      };

      template <typename real_t, typename n_t>
      struct scale_factor
      {
        BOOST_GPU_ENABLED
        real_t operator()(const n_t &n)
        {
          // see section 5.1.3 in Shima et al. 2009
          return real_t((n*(n-1))/2) / (n/2); 
        }
      };

      // assumes _a have higher multiplicities
      template <typename real_t, typename n_t,
        int   n_a, int   n_b,
        int rw2_a, int rw2_b,
        int rd3_a, int rd3_b,
        int  vt_a, int  vt_b,
        typename tup_t
      >
      BOOST_GPU_ENABLED
      void collide(tup_t tpl, const n_t &col_no)
      {
	// multiplicity change (eq. 12 in Shima et al. 2009)
	thrust::get<n_a>(tpl) -= col_no * thrust::get<n_b>(tpl);

	// wet radius change (eq. 13 in Shima et al. 2009)
	thrust::get<rw2_b>(tpl) = pow(
	  col_no * pow(thrust::get<rw2_a>(tpl), real_t(3./2)) + 
	  pow(thrust::get<rw2_b>(tpl), real_t(3./2))
	  ,
	  real_t(2./3)
	);

	// dry radius change (eq. 13 in Shima et al. 2009)
	thrust::get<rd3_b>(tpl) 
	  = col_no *thrust::get<rd3_a>(tpl) + thrust::get<rd3_b>(tpl);

	// invalidating vt
	thrust::get<vt_b>(tpl) = detail::invalid;

	// TODO: kappa, chemistry (only if enabled)
      }

      template <typename real_t, typename n_t>
      struct collider
      {
        // read-only parameters
        typedef thrust::tuple<
          real_t,                       // random number (u01)
          real_t,                       // scaling factor
          thrust_size_t, thrust_size_t, // ix
          thrust_size_t, thrust_size_t, // off (index within cell)
          real_t                       // dv
        > tpl_ro_t;
        enum { u01_ix, scl_ix, ix_a_ix, ix_b_ix, off_a_ix, off_b_ix, dv_ix };

        // read-write parameters = return type
        typedef thrust::tuple<
          n_t,           n_t,           // n   (multiplicity)
          real_t,        real_t,        // rw2 (wet radius squared)
          real_t,        real_t,        // vt  (terminal velocity)
          real_t,        real_t,        // rd3 (dry radius cubed)
          real_t,        real_t         // number of collisions (output); same vector as u01!
        > tpl_rw_t;
        enum { n_a_ix, n_b_ix, rw2_a_ix, rw2_b_ix, vt_a_ix, vt_b_ix, rd3_a_ix, rd3_b_ix, col_a_ix, col_b_ix };

        // read-only parameters passed to the calc function
        typedef thrust::tuple<
          real_t,                      // rhod (dry air density)
          real_t                       // eta (dynamic viscosity)
        > tpl_ro_calc_t;
        enum { rhod_ix, eta_ix };

        const real_t dt;
        const kernel_base<real_t, n_t> *p_kernel;

        //ctor
        collider(const real_t &dt, kernel_base<real_t, n_t> *p_kernel) : dt(dt), p_kernel(p_kernel) {}

        template <class tup_ro_rw_t>
        BOOST_GPU_ENABLED
        void operator()(tup_ro_rw_t tpl_ro_rw)
        {
          const tpl_ro_t &tpl_ro(thrust::get<0>(tpl_ro_rw));
          const tpl_rw_t &tpl_rw(thrust::get<1>(tpl_ro_rw));
          const tpl_ro_calc_t &tpl_ro_calc(thrust::get<2>(tpl_ro_rw));

          // sanity checks
#if !defined(__NVCC__)
          assert(thrust::get<ix_a_ix>(tpl_ro) + 1 == thrust::get<ix_b_ix>(tpl_ro));
#endif
 
          // checking if valid candidates for collision
          {
            const thrust_size_t &cix_a = thrust::get<ix_a_ix>(tpl_ro) - thrust::get<off_a_ix>(tpl_ro);

            // only every second droplet within a cell
            if (cix_a % 2 != 0) return;

            const thrust_size_t &cix_b = thrust::get<ix_b_ix>(tpl_ro) - thrust::get<off_b_ix>(tpl_ro);

            // only droplets within the same cell
            if (cix_a != cix_b - 1)
            {
              thrust::get<col_a_ix>(thrust::get<1>(tpl_ro_rw)) = real_t(0.);
              return;
            }
          }

          //wrap the tpl_rw and tpl_ro_calc tuples to pass it to kernel
          tpl_calc_wrap<real_t,n_t> tpl_wrap(tpl_rw, tpl_ro_calc);

          // computing the probability of collision
          real_t prob = dt / thrust::get<dv_ix>(tpl_ro)
            * thrust::get<scl_ix>(tpl_ro)
            * p_kernel->calc(tpl_wrap);
  
          n_t col_no = n_t(prob); //number of collisions between the pair; rint?

          // comparing the upscaled probability with a random number and returning if unlucky
          if (thrust::get<u01_ix>(tpl_ro) < prob - col_no) ++col_no;
          if(col_no == 0) 
          {
            thrust::get<col_a_ix>(thrust::get<1>(tpl_ro_rw)) = real_t(0.);
            thrust::get<col_b_ix>(thrust::get<1>(tpl_ro_rw)) = real_t(0.);
            return;
          }

#if !defined(__NVCC__)
          using std::min;
#endif
          // performing the coalescence event if lucky
          // note: >= causes equal-multiplicity collisions to result in flagging for recycling
          if (thrust::get<n_a_ix>(tpl_rw) >= thrust::get<n_b_ix>(tpl_rw)) 
          {
            if(thrust::get<n_b_ix>(tpl_rw) > 0) 
              col_no = min( col_no, n_t(thrust::get<n_a_ix>(tpl_rw) / thrust::get<n_b_ix>(tpl_rw)));
            collide<real_t, n_t,
                n_a_ix,   n_b_ix,
              rw2_a_ix, rw2_b_ix,
              rd3_a_ix, rd3_b_ix,
               vt_a_ix,  vt_b_ix
            >(thrust::get<1>(tpl_ro_rw), col_no);
            thrust::get<col_b_ix>(thrust::get<1>(tpl_ro_rw)) = real_t(na_ge_nb); // col vector for the second in a pair stores info on which one has greater multiplicity
          }
          else
          {
            if(thrust::get<n_a_ix>(tpl_rw) > 0) 
              col_no = min( col_no, n_t(thrust::get<n_b_ix>(tpl_rw) / thrust::get<n_a_ix>(tpl_rw)));
            collide<real_t, n_t,
                n_b_ix,   n_a_ix,
              rw2_b_ix, rw2_a_ix,
              rd3_b_ix, rd3_a_ix,
               vt_b_ix,  vt_a_ix
            >(thrust::get<1>(tpl_ro_rw), col_no);
            thrust::get<col_b_ix>(thrust::get<1>(tpl_ro_rw)) = real_t(nb_gt_na); // col vector for the second in a pair stores info on which one has greater multiplicity
          }
          thrust::get<col_a_ix>(thrust::get<1>(tpl_ro_rw)) = real_t(col_no); // col vector for the first in a pair stores info on number of collisions
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::coal(const real_t &dt)
    {   
      // prerequisites
      hskpng_shuffle_and_sort(); // to get random neighbours by default
      hskpng_count();            // no. of super-droplets per cell 
      
      // placing scale_factors in count_mom (of size count_n!)
      thrust::transform(
        count_num.begin(), count_num.begin() + count_n, // input - 1st arg
        count_mom.begin(),                              // output
        detail::scale_factor<real_t, n_t>()
      );

      // references to tmp data
      thrust_device::vector<real_t> 
        &scl(tmp_device_real_cell), // scale factor for probablility
        &col(tmp_device_real_part); // number of collisions, used in chemistry, NOTE: it's the same as u01, so it overwrites already used random numbers
                                    // 1st one of a pair stores number of collisions, 2nd one stores info on which one has greater multiplicity
      thrust_device::vector<thrust_size_t> 
        &off(tmp_device_size_cell); // offset for getting index of particle within a cell

      // laying out scale factor onto ijk grid
      // fill with 0s if not all cells will be updated in the following copy
      if(count_n!=n_cell)  thrust::fill(scl.begin(), scl.end(), real_t(0.));
      
      thrust::copy(
        count_mom.begin(),                    // input - begin
        count_mom.begin() + count_n,          // input - end
        thrust::make_permutation_iterator(    // output
          scl.begin(),                        // data
          count_ijk.begin()                   // permutation
        )
      );  

      // cumulative sum of count_num -> (i - cumsum(ijk(i))) gives droplet index in a given cell
      // fill with 0s if not all cells will be updated in the following copy
      if(count_n!=n_cell)  thrust::fill(off.begin(), off.end(), real_t(0.));
      thrust::copy(
        count_num.begin(), 
        count_num.begin() + count_n, 
        thrust::make_permutation_iterator(    // output
          off.begin(),                        // data
          count_ijk.begin()                   // permutation
        )
      );
      thrust::exclusive_scan( 
        off.begin(), off.end(),
        off.begin()
      );

      // tossing n_part/2 random numbers for comparing with probability of collisions in a pair of droplets
      rand_u01(n_part);

      // colliding
      typedef thrust::permutation_iterator<
        typename thrust_device::vector<thrust_size_t>::iterator,
        typename thrust_device::vector<thrust_size_t>::iterator
      > pi_size_t;

      typedef thrust::permutation_iterator<
        typename thrust_device::vector<real_t>::iterator,
        typename thrust_device::vector<thrust_size_t>::iterator
      > pi_real_t;

      typedef  typename thrust_device::vector<real_t>::iterator i_real_t;

      typedef thrust::permutation_iterator<
        typename thrust_device::vector<n_t>::iterator,
        typename thrust_device::vector<thrust_size_t>::iterator
      > pi_n_t;

      typedef thrust::zip_iterator<
        thrust::tuple< 
          i_real_t,                                                // u01
          pi_real_t,                                               // scl
          thrust::counting_iterator<thrust_size_t>,                // ix_a
          thrust::counting_iterator<thrust_size_t>,                // ix_b
          pi_size_t, pi_size_t,                                    // off_a & off_b
          pi_real_t                                                // dv
        >
      > zip_ro_t;

      typedef thrust::zip_iterator<
        thrust::tuple< 
          pi_n_t,    pi_n_t,    // n_a,   n_b
          pi_real_t, pi_real_t, // rw2_a, rw2_b
          pi_real_t, pi_real_t, // vt_a,  vt_b
          pi_real_t, pi_real_t, // vt_a,  vt_b
          i_real_t, i_real_t    // col_a, col_b 
        >
      > zip_rw_t;

      typedef thrust::zip_iterator<
        thrust::tuple<
          pi_real_t,  // rhod
          pi_real_t   // eta
        >
      > zip_ro_calc_t;    //read-only parameters passed to the calc() function, later also epsilon and Re_lambda

      zip_ro_t zip_ro_it(
        thrust::make_tuple(
          // u01
          u01.begin(),
          // scl
          thrust::make_permutation_iterator(scl.begin(), sorted_ijk.begin()), 
          // ix
          zero,
          zero+1,
          // cid
          thrust::make_permutation_iterator(off.begin(), sorted_ijk.begin()), 
          thrust::make_permutation_iterator(off.begin(), sorted_ijk.begin())+1,
          // dv
          thrust::make_permutation_iterator(dv.begin(), sorted_ijk.begin())
        )
      );

      zip_ro_calc_t zip_ro_calc_it(
        thrust::make_tuple(
          // rhod
          thrust::make_permutation_iterator(rhod.begin(), sorted_ijk.begin()),
          // eta
          thrust::make_permutation_iterator(eta.begin(), sorted_ijk.begin())
        )
      );

      zip_rw_t zip_rw_it(
        thrust::make_tuple(
          // multiplicity
          thrust::make_permutation_iterator(n.begin(),   sorted_id.begin()),  
          thrust::make_permutation_iterator(n.begin(),   sorted_id.begin())+1,
          // wet radius squared
          thrust::make_permutation_iterator(rw2.begin(), sorted_id.begin()), 
          thrust::make_permutation_iterator(rw2.begin(), sorted_id.begin())+1,  
          // terminal velocity
          thrust::make_permutation_iterator(vt.begin(),  sorted_id.begin()), 
          thrust::make_permutation_iterator(vt.begin(),  sorted_id.begin())+1,  
          // dry radius cubed
          thrust::make_permutation_iterator(rd3.begin(), sorted_id.begin()), 
          thrust::make_permutation_iterator(rd3.begin(), sorted_id.begin())+1,
          // output
          col.begin(),  // number of collisions
          col.begin()+1 // which one has greater multiplicity
        )
      );

      thrust::for_each(
        thrust::make_zip_iterator(thrust::make_tuple(zip_ro_it, zip_rw_it, zip_ro_calc_it)),
        thrust::make_zip_iterator(thrust::make_tuple(zip_ro_it, zip_rw_it, zip_ro_calc_it)) + n_part - 1,
        detail::collider<real_t, n_t>(dt, p_kernel)
      );

      // add masses of chemicals
      if(opts_init.chem_switch)
      {
        for(int i=0; i<chem_all; ++i)
          thrust::for_each(
            thrust::make_zip_iterator(thrust::make_tuple(
              thrust::make_permutation_iterator(chem_bgn[i], sorted_id.begin()),
              thrust::make_permutation_iterator(chem_bgn[i], sorted_id.begin())+1,
              col.begin(),
              col.begin()+1
            )),
            thrust::make_zip_iterator(thrust::make_tuple(
              thrust::make_permutation_iterator(chem_bgn[i], sorted_id.begin()),
              thrust::make_permutation_iterator(chem_bgn[i], sorted_id.begin())+1,
              col.begin(),
              col.begin()+1
            )) + n_part -1,
            detail::summator()
          );
      }
      // add kappas
      if(opts_init.dry_distros.size() > 1)
      {
        thrust::for_each(
          thrust::make_zip_iterator(thrust::make_tuple(
            thrust::make_permutation_iterator(kpa.begin(), sorted_id.begin()),    // values to be added
            thrust::make_permutation_iterator(kpa.begin(), sorted_id.begin())+1,
            thrust::make_permutation_iterator(rd3.begin(), sorted_id.begin()),    // weighting factors
            thrust::make_permutation_iterator(rd3.begin(), sorted_id.begin())+1,
            col.begin(),                                                          // collision information
            col.begin()+1
          )),
          thrust::make_zip_iterator(thrust::make_tuple(
            thrust::make_permutation_iterator(kpa.begin(), sorted_id.begin()),
            thrust::make_permutation_iterator(kpa.begin(), sorted_id.begin())+1,
            thrust::make_permutation_iterator(rd3.begin(), sorted_id.begin()),
            thrust::make_permutation_iterator(rd3.begin(), sorted_id.begin())+1,
            col.begin(),
            col.begin()+1
          )) + n_part -1,
          detail::weighted_summator<real_t>()
        );
      }
    }
  };  
};
