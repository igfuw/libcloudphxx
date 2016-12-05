#pragma once

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {   
      // evaluate a unary function and multiply it
      template <typename real_t>
      struct eval_and_mul
      {   
        const common::unary_function<real_t> &fun;
        const real_t &mul;

        // ctor
        eval_and_mul(
          const common::unary_function<real_t> &fun, 
          const real_t &mul
        )
          : fun(fun), mul(mul)
        {}

        BOOST_GPU_ENABLED
        real_t operator()(real_t x) const
        {
            return mul * fun(x);
        }
      }; 

      // evaluate a unary function and add to it
      template <typename real_t>
      struct eval_and_add
      {   
        const common::unary_function<real_t> &fun;
        const real_t &add;

        // ctor
        eval_and_add(
          const common::unary_function<real_t> &fun, 
          const real_t &add
        )
          : fun(fun), add(add)
        {}

        BOOST_GPU_ENABLED
        real_t operator()(real_t x) const
        {
            return add + fun(x);
        }
      }; 
 
      /// @brief returns real_t(exp(3*x))
      template <typename real_t>
      struct exp3x
      {   
        BOOST_GPU_ENABLED 
        real_t operator()(real_t x)  
        {
          return exp(3*x); 
        }
      }; 
    };
  };
};


