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
      // workaround for a clang bug
      typedef std::pair<const real_t*, const ptrdiff_t*> arrinfo_t;

      virtual void init() { assert(false); }  
      virtual void step() { assert(false); }  
      virtual void sync_e2l(
        arrinfo_t rhod_th,
        arrinfo_t rhod_rv,
        arrinfo_t rhod = arrinfo_t(NULL, NULL)
      ) { assert(false); }
      virtual void sync_l2e(
        arrinfo_t rhod_th,
        arrinfo_t rhod_rv
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

      //
      // pair.first -> blitz::Array.dataZero()
      // pair.second -> blitz::Array.stride().data()
      void sync_e2l(
        typename parent_t::arrinfo_t rhod_th,
        typename parent_t::arrinfo_t rhod_rv,
        typename parent_t::arrinfo_t rhod = typename parent_t::arrinfo_t(NULL, NULL)
      );
      void sync_l2e(
        typename parent_t::arrinfo_t rhod_th,
        typename parent_t::arrinfo_t rhod_rv
      );

      void init(); // TODO: explain why init not within constructor? (e.g. NaNs in Eulerian part?)
      void step();
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
