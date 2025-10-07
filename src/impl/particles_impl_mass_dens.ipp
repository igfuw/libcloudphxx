namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {   
      template <typename real_t>
      struct mass_dens_estimator //: thrust::unary_function<const thrust::tuple<real_t, real_t, real_t>&, real_t>
      {

        real_t rad, sig0, xp;

        mass_dens_estimator(real_t rad, real_t sig0, real_t xp) : rad(rad), sig0(sig0), xp(xp) {}

        BOOST_GPU_ENABLED
        real_t operator()(const thrust::tuple<real_t, real_t, real_t> &tpl)
        {
#if !defined(__NVCC__)
          using std::pow;
          using std::exp;
#endif
          const real_t n = thrust::get<0>(tpl);
          const real_t x = thrust::get<1>(tpl);
          const real_t sig = sig0 / pow(thrust::get<2>(tpl), real_t(0.2));
          return n / sig * pow(x, 3 * xp) * exp ( - pow ( (log (pow(x, xp)) - log (rad)) / sig, 2) / 2.); // for wet xp = 1/2, for dry xp = 1/3
        }
      };
    };

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::mass_dens_estim(
      const typename thrust_device::vector<real_t>::iterator &vec_bgn,
      const real_t radius, const real_t sigma0, const real_t power
    )
    {
      assert(selected_before_counting); //TODO: force moms_all() before mass density estimation?

      // same as above
      auto n_filtered_g = tmp_device_real_part.get_guard();
      thrust_device::vector<real_t> &n_filtered = n_filtered_g.get();

      // number of SD in each cell casted to real_t
      thrust_device::vector<real_t> &count_num_real_t(tmp_device_real_cell);

      // get number of SD in each cell
      hskpng_count();

      // cast it to real_t
      thrust::copy(
        count_num.begin(),
        count_num.end(),
        count_num_real_t.begin()
      );


      typedef thrust::permutation_iterator<
        typename thrust_device::vector<real_t>::const_iterator,
        typename thrust_device::vector<thrust_size_t>::iterator
      > pi_t;

      typedef thrust::zip_iterator<thrust::tuple<pi_t, pi_t, pi_t> > zip_it_t;

      thrust::pair<
        thrust_device::vector<thrust_size_t>::iterator,
        typename thrust_device::vector<real_t>::iterator
      > it_pair = thrust::reduce_by_key(
        // input - keys
        sorted_ijk.begin(), sorted_ijk.end(),
        // input - values
        thrust::make_transform_iterator(
          zip_it_t(thrust::make_tuple(
            pi_t(n_filtered.begin(),          sorted_id.begin()),
            pi_t(vec_bgn,                     sorted_id.begin()),
            pi_t(count_num_real_t.begin(),    sorted_ijk.begin())
          )),
          detail::mass_dens_estimator<real_t>(radius, sigma0, power)
        ),
        // output - keys
        count_ijk.begin(),
        // output - values
        count_mom.begin()
      );

      count_n = it_pair.first - count_ijk.begin();
      assert(count_n > 0 && count_n <= n_cell);

#if !defined(__NVCC__)
      using std::sqrt;
#endif

      real_t prefactor = 4. / 3. * ( common::moist_air::rho_w<real_t>() / si::kilograms * si::cubic_metres ) *
#if !defined(__NVCC__)
        sqrt(pi<real_t>() / 2.);
#else
        sqrt(CUDART_PI / 2.);
#endif
      namespace arg = thrust::placeholders;

      //multiply by prefactor and divide by dv
      thrust::transform(
        count_mom.begin(), count_mom.begin() + count_n,     // input - first arg
        thrust::make_permutation_iterator(                  // input - second arg
          dv.begin(),
          count_ijk.begin()
        ),
        count_mom.begin(),                                  // output (in place)
        prefactor * arg::_1 / arg::_2
      );
    }
  };
};
