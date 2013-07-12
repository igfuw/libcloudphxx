#pragma once 

#include <cassert>
#include <memory>

#include "../common/unary_function.hpp"

namespace libcloudphxx
{
  namespace common
  {
    namespace prtcls
    {
      // to make inclusion of Thrust not neccesarry
      enum {cpp, omp, cuda};

      // to allow storing instances for multiple backends in one container/pointer
      template <typename real_t>
      class particles_proto // TODO: rename to any?
      {
        public: 
        virtual void init(unary_function<real_t> *pdf) { assert(false); }  
      };  

      // prototype of what's implemented in the .tpp file
      template <typename real_t, int thrust_device_system>
      class particles : public particles_proto<real_t>
      {
        // pimpl stuff
        struct impl;
        std::auto_ptr<impl> pimpl;
    
        // the public API
        public:  
        particles(real_t sd_conc_mean, int nx, int ny, int nz); // ctor
        void init(unary_function<real_t> *pdf);
      };

      // to be explicitely instantiated
      template <typename real_t>
      particles_proto<real_t> *factory(
        int backend,
        real_t sd_conc_mean, 
        int nx = 0, int ny = 0, int nz = 0
      ); 
    };
  };
};
