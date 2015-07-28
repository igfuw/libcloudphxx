import sys
sys.path.insert(0, "../../bindings/python/")

from libcloudphxx import lgrngn

from numpy import array as arr_t, frombuffer, repeat, zeros, float64

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
opts_init.kernel = lgrngn.kernel_t.geometric
opts_init.dt = 1
opts_init.sd_conc = 64
opts_init.n_sd_max = 512
opts_init.rng_seed = 396

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

print "chem_switch = ", opts_init.chem_switch
print "dt =", opts_init.dt
print "sstp_cond =", opts_init.sstp_cond
print "sstp_coal =", opts_init.sstp_coal
print "sstp_chem =", opts_init.sstp_chem 

print "kernel =", opts_init.kernel 
print "sd_conc =", opts_init.sd_conc
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

prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)
prtcls.step_sync(opts, th, rv, rhod)
rain = prtcls.step_async(opts)
prtcls.diag_dry_rng(0.,1.)
prtcls.diag_wet_rng(0.,1.)
prtcls.diag_dry_mom(1)
prtcls.diag_wet_mom(1)
prtcls.diag_all()
#prtcls.diag_chem(lgrngn.chem_species_t.OH)
prtcls.diag_sd_conc()
assert frombuffer(prtcls.outbuf()) == opts_init.sd_conc # parcel set-up

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
prtcls.init(th, rv, rhod) #TODO: test passing rhoCx, rhoCy, rhoCz here 
prtcls.diag_sd_conc()
assert len(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.nx
assert (frombuffer(prtcls.outbuf()) > 0).all()
assert sum(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.nx * opts_init.sd_conc

# 3arg variant (const rho)
prtcls.step_sync(opts, th, rv) #TODO: this should fail as no Courants were passed in init!
rain = prtcls.step_async(opts)

# 4arg variant (var rho - the one used in 0D)
prtcls.step_sync(opts, th, rv, rhod) #TODO: this should fail as no Courants were passed in init!
rain = prtcls.step_async(opts)

#TODO: test profile vs. 2D array

# 3D
rhod = arr_t([rhod, rhod])
th   = arr_t([th,   th  ])
rv   = arr_t([rv,   rv  ])

opts_init.ny = 2
opts_init.dy = 10
opts_init.y1 = opts_init.ny * opts_init.dy

prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod) #TODO: test passing rhoCx, rhoCy, rhoCz here
prtcls.diag_sd_conc()
assert len(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.nx * opts_init.ny
assert (frombuffer(prtcls.outbuf()) > 0).all()
assert sum(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.nx * opts_init.ny * opts_init.sd_conc

# 3 arg variant
prtcls.step_sync(opts, th, rv) #TODO: should fail with no Courant in init
rain = prtcls.step_async(opts)

# 4 arg variant
prtcls.step_sync(opts, th, rv, rhod) #TODO: should fail with no Courants in init
rain = prtcls.step_async(opts)

# 6 arg variant
rhoCx = zeros((opts_init.nx+1, opts_init.ny+0, opts_init.nz+0), dtype=float64) #TODO: these dimensions are not checked...
rhoCy = zeros((opts_init.nx+0, opts_init.ny+1, opts_init.nz+0), dtype=float64)
rhoCz = zeros((opts_init.nx+0, opts_init.ny+0, opts_init.nz+1), dtype=float64)
prtcls.step_sync(opts, th, rv, rhoCx, rhoCy, rhoCz)
rain = prtcls.step_async(opts)

#TODO: test profile vs. 3D array
