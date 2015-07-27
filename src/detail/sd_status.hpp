#pragma once 

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      enum sd_stat_t { active, inactive, to_init, to_rcyc}; 

      // status getters
      struct is_active_or_to_init : std::unary_function<sd_stat_t, bool>
      {
        BOOST_GPU_ENABLED
        bool operator()(const sd_stat_t &stat)
        {
          return (stat == active || stat == to_init);
        }
      };

      struct is_active : std::unary_function<sd_stat_t, bool>
      {
        BOOST_GPU_ENABLED
        bool operator()(const sd_stat_t &stat)
        {
          return (stat==active);
        }
      };

      struct is_to_init : std::unary_function<sd_stat_t, bool>
      {
        BOOST_GPU_ENABLED
        bool operator()(const sd_stat_t &stat)
        {
          return (stat == to_init);
        }
      };

      // manipualtors
      template <typename real_t>
      struct deactivate
      {   
        BOOST_GPU_ENABLED
        sd_stat_t operator()(const real_t &)
        {
          return inactive;
        }
      };

      template <typename real_t>
      struct activate
      {   
        BOOST_GPU_ENABLED
        sd_stat_t operator()(const real_t &)
        {
          return active;
        }
      };

      // comparison
      struct active_first : std:binary_function<sd_stat_t, sd_stat_t, bool>
      {
        BOOST_GPU_ENABLED
        bool operator()(const sd_stat_t &a, const sd_stat_t &b) { return is_active(a);}         
      }
    };
  };
};
