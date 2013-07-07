#include <cassert>

// to make inclusion of Thrust not neccesarry
enum {cpp, openmp, cuda};

// to allow storing instances for multiple backends in one container/pointer
template <typename real_t>
class particles_proto
{
  public: 
  virtual void func() { assert(false); }   
};  

// prototype of what's implemented in the .tpp file
template <typename real_t, int thrust_device_system>
class particles : public particles_proto<real_t>
{
  public:  
  void func();
};
