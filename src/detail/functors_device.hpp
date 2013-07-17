// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {

      /// @brief returns real_t(exp(3*x))
      template <typename real_t>
      struct exp3x
      { 
	__device__ 
	real_t operator()(real_t x) 
	{ 
	  return exp(3*x); 
	} 
      };


      /// @brief returns ret_t(x/c) 
      template <typename arg_t, typename ret_t>
      struct divide_by_constant_and_cast
      {
        arg_t c;
        divide_by_constant_and_cast(arg_t c) : c(c) {}

        __device__
        ret_t operator()(arg_t x) 
        { 
          return ret_t(x/c); 
        }
      };


      /// @brief returns i+j*n (ravels {i,j} index pairs into a single index)
// TODO: what about 3D
// TODO: probably can be done using placeholders!
      struct ravel_indices
      {
        int n;
        ravel_indices(int n) : n(n) {}
      
        __device__
        int operator()(int i, int j) 
        { 
          return i + j * n; 
        }
      };
    };
  };
};
