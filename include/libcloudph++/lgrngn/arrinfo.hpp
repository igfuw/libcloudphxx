#pragma once 

#include <libcloudph++/lgrngn/extincl.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    // helper struct to ease passing n-dimensional arrays
//<listing>
    template <typename real_t>
    struct arrinfo_t
    {
      // member fields:
      real_t * const dataZero;
      const ptrdiff_t *strides;

//</listing>
      const std::vector<ptrdiff_t> strvec; // see alt. ctor below

      // ctors
      arrinfo_t()
        : dataZero(NULL), strides(NULL) 
      {} 

      arrinfo_t(real_t * const dataZero, const ptrdiff_t *strides) 
        : dataZero(dataZero), strides(strides) 
      {} 

      // methods
      bool is_null() const { return dataZero==NULL || strides==NULL; }

      // alternative usage with local storage of the strides
      arrinfo_t(real_t * const dataZero, const std::vector<ptrdiff_t> &strvec) :
        strvec(strvec), dataZero(dataZero), strides(&strvec[0])
      {}

      // non-default copy ctor handling both the original and alternative usage
      arrinfo_t(const arrinfo_t &ai) :
        strvec(ai.strvec), 
        dataZero(ai.dataZero), 
        strides(!strvec.empty() ? &strvec[0] : ai.strides)
      {}
    };
  };
};
