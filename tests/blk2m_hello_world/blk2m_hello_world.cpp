#include <memory>
#include <set>

#include <blitz/tv2fastiter.h> // otherwise Clang fails in debug mode 
#include <blitz/array.h>

#include <libcloudph++/blk_2m/options.hpp>

#include <libcloudph++/blk_2m/rhs_cellwise.hpp>
#include <libcloudph++/blk_2m/rhs_columnwise.hpp>

namespace b2m = libcloudphxx::blk_2m;
using real_t = double;
const int 
  n_dims = 2,
  n_halo = 1,
  nx = 4, 
  nz = 4;
const real_t dt = 1, dz = 10;

// MicroHH-style (restrict + ijk Fortran index order)
using ptr_t = std::unique_ptr<real_t>; // not using restrict as Blitz++ anyhow assumes it

blitz::Array<real_t, n_dims> bz(const ptr_t &ptr)
{
  return blitz::Array<real_t, n_dims>(
    ptr.get(), 
    blitz::shape(
      nx + 2 * n_halo, 
      nz + 2 * n_halo
    ), 
    blitz::neverDeleteData,
    blitz::ColumnMajorArray<n_dims>()
  ).reindex({
    -n_halo, 
    -n_halo
  });
}

blitz::Range 
  i_all(0, nx-1), 
  j_all(0, nz-1);

int main()
{
  ptr_t
    rho_d,
    rr, nr, 
    dot_rr, dot_nr
  ;
  
  // memory allocation
  for (auto arr : std::set<ptr_t*>({ &rho_d, &rr, &nr, &dot_rr, &dot_nr }))
    arr->reset(new real_t[ (nx+2*n_halo) * (nz+2*n_halo) ]);

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

  // some meaningfull values
  bz_rho_d(i_all, j_all) = 1;
  bz_rr(i_all, nz-1) = 1e-5;  
  bz_nr(i_all, nz-1) = 1e3;  

  std::cerr << bz_rr << std::endl;

  // calling the blk_2m methods from libcloudph++ 
  for (int i=0; i < nx; ++i)
  {
    // columns of dotted fields to be modified
    auto 
      bz_dot_rr_ij = bz_dot_rr(i, j_all),
      bz_dot_nr_ij = bz_dot_nr(i, j_all);
    b2m::rhs_columnwise(
      opts,
      bz_dot_rr_ij,
      bz_dot_nr_ij,
      bz_rho_d  (i, j_all),
      bz_rr     (i, j_all),
      bz_nr     (i, j_all),
      dt,
      dz
    );
  }

  std::cerr << bz_dot_rr << std::endl;

  //b2m::rhs_cellwise(); // TODO
}
