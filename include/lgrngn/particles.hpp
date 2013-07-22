#pragma once 

#include <cassert>
#include <cstddef> // ptrdiff_t
#include <memory>


#include <boost/ptr_container/ptr_unordered_map.hpp>

#include "../common/unary_function.hpp"

namespace libcloudphxx
{
  namespace lgrngn
  {
    using libcloudphxx::common::unary_function;

    // to make inclusion of Thrust not neccesarry
    enum {cpp, omp, cuda};

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
    };  

    template <typename real_t>
    struct opts_t 
    {   
      typedef boost::ptr_unordered_map<real_t, unary_function<real_t> > dry_distros_t;

      int nx, ny, nz;
      real_t sd_conc_mean; 
      real_t dx, dy, dz;
      dry_distros_t dry_distros;
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
