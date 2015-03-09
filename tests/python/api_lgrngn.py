import sys
sys.path.insert(0, "../../bindings/python/")

from libcloudphxx import lgrngn

from numpy import array as arr_t, frombuffer

from math import exp, log, sqrt, pi

def lognormal(lnr):
  mean_r = .04e-6 / 2
  stdev  = 1.4
  n_tot  = 60e6
  return n_tot * exp(
    -pow((lnr - log(mean_r)), 2) / 2 / pow(log(stdev),2)
  ) / log(stdev) / sqrt(2*pi);

opts_init = lgrngn.opts_init_t()
kappa = .61
opts_init.dry_distros = {kappa:lognormal}

print "nx =", opts_init.nx
print "ny =", opts_init.ny
print "nz =", opts_init.nz

print "dx =", opts_init.dx
print "dy =", opts_init.dy
print "dz =", opts_init.dz

print "x0 =", opts_init.x0
print "y0 =", opts_init.y0
print "z0 =", opts_init.z0

print "x1 =", opts_init.x1
print "y1 =", opts_init.y1
print "z1 =", opts_init.z1

print "dt =", opts_init.dt
print "sstp_cond =", opts_init.sstp_cond
print "sstp_coal =", opts_init.sstp_coal
print "sstp_chem =", opts_init.sstp_chem 

print "kernel =", opts_init.kernel 
assert opts_init.kernel == lgrngn.kernel_t.geometric
print "sd_conc_mean =", opts_init.sd_conc_mean
print "chem_rho =", opts_init.chem_rho

print lgrngn.backend_t.OpenMP
print lgrngn.backend_t.CUDA
backend = lgrngn.backend_t.serial

opts = lgrngn.opts_t()
print "adve =", opts.adve
print "sedi =", opts.sedi
print "cond =", opts.cond
print "coal =", opts.coal
print "chem =", opts.chem
print "RH_max =", opts.RH_max

opts.chem_gas = {
  lgrngn.chem_species_t.SO2  : 44,
  lgrngn.chem_species_t.O3   : 44,
  lgrngn.chem_species_t.H2O2 : 44
}

print "chem_gas[SO2] = ", opts.chem_gas[lgrngn.chem_species_t.SO2]
print "chem_gas = ", opts.chem_gas

# 0D (parcel)
rhod = arr_t([  1.])
th   = arr_t([300.])
rv   = arr_t([  0.01])

print rhod.dtype, arr_t([  1.])

prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)
prtcls.step_sync(opts, th, rv, rhod)
rain = prtcls.step_async(opts)
prtcls.diag_dry_rng(0.,1.)
prtcls.diag_wet_rng(0.,1.)
prtcls.diag_dry_mom(1)
prtcls.diag_wet_mom(1)
prtcls.diag_all()
prtcls.diag_chem(lgrngn.chem_species_t.OH)
prtcls.diag_sd_conc()
assert frombuffer(prtcls.outbuf()) == opts_init.sd_conc_mean # parcel set-up

# 1D
rhod = arr_t([  1.,   1.,   1.])
th   = arr_t([300., 300., 300.])
rv   = arr_t([   .01,  .01,  .01])

opts_init.nz = 3 
opts_init.dz = 10
opts_init.z1 = opts_init.nz * opts_init.dz
# TODO... still not implemented in the library itself!
#prtcls = lgrngn.factory(backend, opts_init)
#prtcls.init(th, rv, rhod)
#prtcls.diag_sd_conc()
#assert len(frombuffer(prtcls.outbuf())) == opts_init.nz
#assert (frombuffer(prtcls.outbuf()) > 0).all()
#assert sum(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.sd_conc_mean

# 2D
rhod = arr_t([[  1.,    1.  ],[   1.,     1.  ]])
th   = arr_t([[300.,  300.  ],[ 300.,   300.  ]])
rv   = arr_t([[   .01,   .01],[    .01,    .01]])

opts_init.nz = 2
opts_init.nx = 2
opts_init.dz = 10
opts_init.dx = 10
opts_init.z1 = opts_init.nz * opts_init.dz
opts_init.x1 = opts_init.nx * opts_init.dx

prtcls = lgrngn.factory(backend, opts_init)
assert opts_init.nx == prtcls.opts_init.nx
prtcls.init(th, rv, rhod)
prtcls.diag_sd_conc()
assert len(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.nx
assert (frombuffer(prtcls.outbuf()) > 0).all()
assert sum(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.nx * opts_init.sd_conc_mean

# 3D
# TODO...
