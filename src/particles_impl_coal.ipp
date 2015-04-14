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
      template <typename tup_t, typename real_t, typename n_t,
        int   n_a, int   n_b,
        int rw2_a, int rw2_b,
        int rd3_a, int rd3_b,
        int  vt_a, int  vt_b
      >
      BOOST_GPU_ENABLED
      void collide(tup_t &tpl, const n_t &col_no)
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
          real_t,        real_t         // rd3 (dry radius cubed)
        > tpl_rw_t;
        enum { n_a_ix, n_b_ix, rw2_a_ix, rw2_b_ix, vt_a_ix, vt_b_ix, rd3_a_ix, rd3_b_ix };

        const real_t dt;
        const kernel_base<real_t, n_t> *p_kernel;

        //ctor
        collider(const real_t &dt, kernel_base<real_t, n_t> *p_kernel) : dt(dt), p_kernel(p_kernel) {}

        BOOST_GPU_ENABLED
        tpl_rw_t operator()(const tpl_ro_t &tpl_ro, tpl_rw_t tpl_rw)
        {
          // sanity checks
#if !defined(__NVCC__)
          assert(thrust::get<ix_a_ix>(tpl_ro) + 1 == thrust::get<ix_b_ix>(tpl_ro));
#endif
 
          // checking if valid candidates for collision
          {
            const thrust_size_t &cix_a = thrust::get<ix_a_ix>(tpl_ro) - thrust::get<off_a_ix>(tpl_ro);

            // only every second droplet within a cell
            if (cix_a % 2 != 0) return tpl_rw;

            const thrust_size_t &cix_b = thrust::get<ix_b_ix>(tpl_ro) - thrust::get<off_b_ix>(tpl_ro);

            // only droplets within the same cell
            if (cix_a != cix_b - 1) return tpl_rw;
          }

          //wrap the tpl_rw tuple to pass it to kernel
          tpl_rw_t_wrap<real_t,n_t> tpl_wrap(tpl_rw);

          // computing the probability of collision
          real_t prob = dt / thrust::get<dv_ix>(tpl_ro)
            * thrust::get<scl_ix>(tpl_ro)
            * p_kernel->calc(tpl_wrap);
  
          // sanity check for random sampling validity
//          assert(prob < 1); // TODO: there is a workaround proposed in Shima et al. 2009
          n_t col_no = n_t(prob); //number of collisions between the pair; rint?

          // comparing the upscaled probability with a random number and returning if unlucky
          if (thrust::get<u01_ix>(tpl_ro) < prob - col_no) ++col_no;
          if(col_no == 0) return tpl_rw;

#if !defined(__NVCC__)
          using std::min;
#endif
          // performing the coalescence event if lucky
          // note: >= causes equal-multiplicity collisions to result in flagging for recycling
          if (thrust::get<n_a_ix>(tpl_rw) >= thrust::get<n_b_ix>(tpl_rw)) 
          {
            col_no = min( col_no, n_t(thrust::get<n_a_ix>(tpl_rw) / thrust::get<n_b_ix>(tpl_rw)));
            collide<tpl_rw_t, real_t, n_t,
                n_a_ix,   n_b_ix,
              rw2_a_ix, rw2_b_ix,
              rd3_a_ix, rd3_b_ix,
               vt_a_ix,  vt_b_ix
            >(tpl_rw, col_no);
          }
          else
          {
            col_no = min( col_no, n_t(thrust::get<n_b_ix>(tpl_rw) / thrust::get<n_a_ix>(tpl_rw)));
            collide<tpl_rw_t, real_t, n_t,
                n_b_ix,   n_a_ix,
              rw2_b_ix, rw2_a_ix,
              rd3_b_ix, rd3_a_ix,
               vt_b_ix,  vt_a_ix
            >(tpl_rw, col_no);
          }
          return tpl_rw;          
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::coal(const real_t &dt)
    {   
      // prerequisites
      hskpng_shuffle_and_sort(); // to get random neighbours by default
      hskpng_count();            // no. of particles per cell 
      
      // placing scale_factors in count_mom (of size count_n!)
      thrust::transform(
        count_num.begin(), count_num.begin() + count_n, // input - 1st arg
        count_mom.begin(),                              // output
        detail::scale_factor<real_t, n_t>()
      );

      // references to tmp data
      thrust_device::vector<real_t> 
        &scl(tmp_device_real_cell); // scale factor for probablility
      thrust_device::vector<thrust_size_t> 
        &off(tmp_device_size_cell); // offset for getting index of particle within a cell

      // laying out scale factor onto ijk grid
      thrust::fill(scl.begin(), scl.end(), real_t(0));
      thrust::copy(
        count_mom.begin(),                    // input - begin
        count_mom.begin() + count_n,          // input - end
        thrust::make_permutation_iterator(    // output
          scl.begin(),                        // data
          count_ijk.begin()                   // permutation
        )
      );  

      // cumulative sum of count_num -> (i - cumsum(ijk(i))) gives droplet index in a given cell
      thrust::fill(off.begin(), off.end(), real_t(0));
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
      rand_u01(n_part); // TODO: n_part/2 is enough but how to do it with the logic below???

      // colliding
      typedef thrust::permutation_iterator<
        typename thrust_device::vector<thrust_size_t>::iterator,
        typename thrust_device::vector<thrust_size_t>::iterator
      > pi_size_t;

      typedef thrust::permutation_iterator<
        typename thrust_device::vector<real_t>::iterator,
        typename thrust_device::vector<thrust_size_t>::iterator
      > pi_real_t;

      typedef thrust::permutation_iterator<
        typename thrust_device::vector<n_t>::iterator,
        typename thrust_device::vector<thrust_size_t>::iterator
      > pi_n_t;

      typedef thrust::zip_iterator<
        thrust::tuple< 
          typename thrust_device::vector<real_t>::iterator,        // u01
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
          pi_real_t, pi_real_t  // rd3_a, rd3_b 
        >
      > zip_rw_t;

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
          thrust::make_permutation_iterator(rd3.begin(), sorted_id.begin())+1
        )
      );

      thrust::transform(
        zip_ro_it, zip_ro_it + n_part - 1,  // input - 1st arg
        zip_rw_it,                          // input - 2nd arg
        zip_rw_it,                          // output (in place)
        detail::collider<real_t, n_t>(dt, p_kernel)
      );
    }
  };  
};
