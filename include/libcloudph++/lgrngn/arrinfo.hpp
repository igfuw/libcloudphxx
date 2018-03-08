#pragma once 

#include <libcloudph++/lgrngn/extincl.hpp>

namespace libcloudphxx
{
  namespace lgrngn
  {
    // helper struct to ease passing n-dimensional arrays
    template <typename real_t>
    struct arrinfo_t
    {
      // member fields:
      const std::vector<ptrdiff_t> strvec; // see alt. ctor below
      real_t * const data;
      const ptrdiff_t *strides;


      // ctors
      arrinfo_t()
        : data(NULL), strides(NULL) 
      {} 

      arrinfo_t(real_t * const data, const ptrdiff_t *strides) 
        : data(data), strides(strides) 
      {} 

      // methods
      bool is_null() const { return data==NULL || strides==NULL; }

      // alternative usage with local storage of the strides
      arrinfo_t(real_t * const data, const std::vector<ptrdiff_t> &_strvec) :
        strvec(_strvec), data(data), strides(&strvec[0])
      {}

      // non-default copy ctor handling both the original and alternative usage
      arrinfo_t(const arrinfo_t &ai) :
        strvec(ai.strvec), 
        data(ai.data), 
        strides(!strvec.empty() ? &strvec[0] : ai.strides)
      {}
    };
  };
};
