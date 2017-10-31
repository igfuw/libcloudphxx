namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {   
      // calculate numerical integral with the trapezoidal rule
      // TODO: use thrust
      template<typename real_t>
      real_t integrate(const common::unary_function<real_t> &fun, const real_t &min, const real_t &max, const real_t &bin_size)
      {   
        const int n = (max - min) / bin_size; //no of bins
        real_t integral = (fun(min) + fun(max)) / 2.; 

        for(int i=1; i<n; ++i)
          integral += fun(min + i * bin_size); 

        return integral * bin_size;
      }   

      // calculate cumulative distribution function
      template<typename real_t, typename vec_t>
      void calc_CDF(const common::unary_function<real_t> &fun, const real_t &min, const real_t &max, const real_t &bin_size, vec_t &vec)
      {   
        const int n = (max - min) / bin_size + 1; //no of points at which cdf will be calculated
        vec.resize(n);

        namespace arg = thrust::placeholders;
        // fill vec with fun values at n points
        thrust::transform(
          thrust::make_transform_iterator(thrust::make_counting_iterator(0), min + bin_size * arg::_1),
          thrust::make_transform_iterator(thrust::make_counting_iterator(n), min + bin_size * arg::_1),
          vec.begin(), eval_and_oper<real_t>(fun, 1));

        // calculate CDF
        thrust::inclusive_scan(vec.begin(), vec.end(), vec.begin());
    
        // normalize CDF     
        thrust::transform(vec.begin(), vec.end(), vec.begin(), arg::_1 / vec.back());
      }   
    };  

    template <typename real_t, backend_t device>
    void particles_t<real_t, device>::impl::dist_analysis(const common::unary_function<real_t> *n_of_lnrd_stp)
    {
      const real_t rd_min_init = 1e-11, rd_max_init = 1e-3;
      rd_min = rd_min_init;
      rd_max = rd_max_init;

      boost::uintmax_t max_iter = 1e6; // max number of iterations when searching for maxima/roots

      // analysing initial distribution
      // has to have only single maximum, no local minima
      // TODO: extend this to more robust initial function!!!
      std::pair<real_t, real_t> init_distr_max; // [ln(position of distribution's maximum), -function value at maximum]
      init_distr_max = boost::math::tools::brent_find_minima(detail::eval_and_oper<real_t>(*n_of_lnrd_stp, -1), log(rd_min), log(rd_max), 32, max_iter);

      real_t init_dist_bound_value = -init_distr_max.second / 1e4; // value of the distribution at which we bind it

      rd_min =
        common::detail::bisect<libcloudphxx::lgrngn::detail::eval_and_oper<real_t>, real_t>(
          detail::eval_and_oper<real_t>(*n_of_lnrd_stp, -init_dist_bound_value, 1),
          real_t(log(rd_min_init)), init_distr_max.first,
          init_dist_prec
        );

      rd_max =
        common::detail::bisect<libcloudphxx::lgrngn::detail::eval_and_oper<real_t>, real_t>(
          detail::eval_and_oper<real_t>(*n_of_lnrd_stp, -init_dist_bound_value, 1),
          init_distr_max.first, real_t(log(rd_max_init)),
          init_dist_prec
        );

      integral = detail::integrate(*n_of_lnrd_stp, rd_min, rd_max, init_dist_prec);

      // calculate cumulative distribution function
      thrust::host_vector<real_t> cdf;

      detail::calc_CDF(*n_of_lnrd_stp, rd_min, rd_max, init_dist_prec, cdf);
    }
  };
};
