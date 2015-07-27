#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      enum sd_stat_t { active, inactive, to_init, to_rcyc}; 

      struct active_or_to_init : std::unary_function<bool, sd_stat_t>
      {
        BOOST_GPU_ENABLED
        bool operator()(const sd_stat_t &stat)
        {
          return (stat == active || stat == to_init);
        }
      };

      struct active : std::unary_function<bool, sd_stat_t>
      {
        BOOST_GPU_ENABLED
        bool operator()(const sd_stat_t &stat)
        {
          return (stat == active);
        }
      };

      struct to_init : std::unary_function<bool, sd_stat_t>
      {
        BOOST_GPU_ENABLED
        bool operator()(const sd_stat_t &stat)
        {
          return (stat == to_init);
        }
      };

      template <typename real_t>
      struct deactivate
      {   
        BOOST_GPU_ENABLED
        detail::sd_stat_t operator()(const real_t &)
        {
          return detail::inactive;
        }
      };

      template <typename real_t>
      struct activate
      {   
        BOOST_GPU_ENABLED
        detail::sd_stat_t operator()(const real_t &)
        {
          return detail::active;
        }
      };
    };
  };
};
