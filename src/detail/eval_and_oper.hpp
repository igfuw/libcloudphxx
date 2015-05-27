#pragma once

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {   
      // evaluate a unary function and multiply it/add to it
      template <typename real_t>
      struct eval_and_oper
      {   
        const common::unary_function<real_t> &fun;
        const real_t &mul;
        const bool operation; // 0 - multiply, 1 - add

        // ctor
        eval_and_oper(
          const common::unary_function<real_t> &fun, 
          const real_t &mul,
          bool operation = 0 
        )
          : fun(fun), mul(mul), operation(operation)
        {}

        real_t operator()(real_t x)  
        {
          if(operation == 0)
            return mul * fun(x);
          else
            return mul + fun(x);
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


