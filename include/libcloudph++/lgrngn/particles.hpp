#pragma once 

#include <cassert>
#include <cstddef> // ptrdiff_t
#include <memory>

#include "options.hpp"

namespace libcloudphxx
{
  namespace lgrngn
  {
    // to make inclusion of Thrust not neccesarry
    enum {cpp, omp, cuda}; // TODO: move to opts?

    // to allow storing instances for multiple backends in one container/pointer
    template <typename real_t>
    struct particles_proto // TODO: rename to any?
    {
      virtual void init(
        const ptrdiff_t* strides,
        real_t *rhod_th,
        real_t *rhod_rv,
        real_t *rhod 
      ) { assert(false); }  

      virtual void step(
        real_t *rhod_th,
        real_t *rhod_rv,
        real_t *rhod = NULL
      ) { assert(false); }  

      virtual void diag(
      ) { assert(false); }

      virtual real_t *outbuf() 
      { 
        assert(false);
        return NULL;
      }
    };  

    // prototype of what's implemented in the .tpp file
    template <typename real_t, int thrust_device_system>
    class particles : public particles_proto<real_t>
    {
      typedef particles_proto<real_t> parent_t;

      // pimpl stuff
      struct impl;
      std::auto_ptr<impl> pimpl;
  
      // the public API
      public:  
      particles(const opts_t<real_t>); // ctor

      // init separated from the ctor as not all data might be available
      void init(
        const ptrdiff_t* strides,
        real_t *rhod_th,
        real_t *rhod_rv,
        real_t *rhod 
      );

      void step(
        real_t *rhod_th,
        real_t *rhod_rv,
        real_t *rhod = NULL
      );

      void diag();

      real_t *outbuf();
    };

    // to be explicitely instantiated
    template <typename real_t>
    struct factory
    {
      static particles_proto<real_t> *make(const int backend, 
	const real_t sd_conc_mean,
	typename opts_t<real_t>::dry_distros_t dry_distros,
	const int = -1, const real_t = 0, 
	const int = -1, const real_t = 0, 
	const int = -1, const real_t = 0
      );
    };
  };
};
