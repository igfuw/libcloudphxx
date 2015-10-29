// vim:filetype=cpp
/** @file
  * @copyright University of Warsaw
  * @section LICENSE
  * GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
  */

#pragma once

// TODO: place contents of ths function into files which use them (apparently, there's no reuse)

#include <libcloudph++/common/unary_function.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    namespace detail
    {
      template <typename real_t>
      struct c_arr_get
      {   
        const real_t * const c_arr; // member field
        c_arr_get(const real_t * const c_arr) : c_arr(c_arr) {} // ctor

        // op invoked by transform
        real_t operator()(const thrust_size_t ix) 
        { 
          return c_arr[ix]; 
        }
      }; 

      template <typename real_t>
      struct c_arr_set
      {   
        real_t *c_arr; // member field
        c_arr_set(real_t *c_arr) : c_arr(c_arr) {} // ctor

        // op invoked by transform
        bool operator()(const thrust_size_t ix, const real_t &val) 
        { 
          c_arr[ix] = val; 
          return true;
        } 
      }; 
    };
  };
};
