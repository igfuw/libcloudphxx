#include <cassert>

enum {cpp, openmp, cuda};

template <typename real_t>
class particles_proto
{
  public: 
  virtual void func()   
  {   
    assert(false);
  }   
};  

template <typename real_t, int thrust_device_system>
class particles : public particles_proto<real_t>
{
  public:  
  void func();
};
