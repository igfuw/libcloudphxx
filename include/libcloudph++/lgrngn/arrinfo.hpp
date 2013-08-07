#pragma once 

#include <cstddef> // ptrdiff_t

namespace libcloudphxx
{
  namespace lgrngn
  {
    // helper struct to ease passing n-dimensional arrays
    template <typename real_t>
    struct arrinfo_t
    {
      // member fields
      real_t * const dataZero;
      const ptrdiff_t *strides;

      // ctors
      arrinfo_t()
        : dataZero(NULL), strides(NULL) 
      {} 

      arrinfo_t(real_t * const dataZero, const ptrdiff_t *strides) 
        : dataZero(dataZero), strides(strides) 
      {} 

      // methods
      bool is_null() const { return dataZero==NULL || strides==NULL; }
    };
  };
};
