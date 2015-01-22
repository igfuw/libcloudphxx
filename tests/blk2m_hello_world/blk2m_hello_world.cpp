#include <memory>
#include <set>

#include <blitz/array.h>

#include <libcloudph++/blk_2m/options.hpp>

#include <libcloudph++/blk_2m/rhs_cellwise.hpp>
#include <libcloudph++/blk_2m/rhs_columnwise.hpp>

namespace b2m = libcloudphxx::blk_2m;
using real_t = double;

const size_t nx = 64, nz = 64;
const real_t dt = 1, dz = 10;

// MicroHH-style (reastrict + ijk Fortran index order)
using ptr_t = std::unique_ptr<real_t restrict>;

blitz::Array<real_t, 2> bz(const ptr_t &ptr)
{
  return blitz::Array<real_t, 2>(
    *ptr, 
    blitz::shape(nx, nz), 
    blitz::FortranArray<2>()
);
}

int main()
{
  ptr_t
    rho_d,
    rr, nr, 
    dot_rr, dot_nr
  ;
  
  // memory allocation
  for (auto arr : std::set<ptr_t*>({ &rho_d, &rr, &nr, &dot_rr, &dot_nr }))
    arr->reset(new real_t[nx*nz]);

  // instantiating the options container
  b2m::opts_t<real_t> opts;  
  
  // using Blitz++ as a way to access pre-allocated data
  auto 
    bz_rho_d  = bz(rho_d),
    bz_rr     = bz(rr),
    bz_nr     = bz(nr),
    bz_dot_rr = bz(dot_rr),
    bz_dot_nr = bz(dot_nr)
  ;

  // calling the blk_2m methods from libcloudph++ 
  b2m::rhs_columnwise(
    opts,
    bz_dot_rr,
    bz_dot_nr,
    bz_rho_d,
    bz_rr,
    bz_nr,
    dt,
    dz
  );

  //b2m::rhs_cellwise(); // TODO
}
