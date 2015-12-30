import sys
sys.path.insert(0, "../../bindings/python/")

from libcloudphxx import lgrngn

from numpy import array as arr_t, frombuffer, repeat, zeros, float64, ones

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
opts_init.terminal_velocity = lgrngn.vt_t.beard76
opts_init.dt = 1
opts_init.sd_conc = 64
opts_init.n_sd_max = 512
opts_init.rng_seed = 396
opts_init.src_dry_distros = {kappa:lognormal}
opts_init.src_sd_conc = 64
opts_init.src_z1 = opts_init.dz

print "nx = ", opts_init.nx
print "ny = ", opts_init.ny
print "nz = ", opts_init.nz
print "dx = ", opts_init.dx
print "dy = ", opts_init.dy
print "dz = ", opts_init.dz
print "x0 = ", opts_init.x0
print "y0 = ", opts_init.y0
print "z0 = ", opts_init.z0
print "x1 = ", opts_init.x1
print "y1 = ", opts_init.y1
print "z1 = ", opts_init.z1

print "sd_conc = ", opts_init.sd_conc
print "RH_max = ", opts_init.RH_max
print "rng_seed = ", opts_init.rng_seed
print "kernel =", opts_init.kernel 
print "sd_conc =", opts_init.sd_conc
print "terminal_velocity =", opts_init.terminal_velocity
print "chem_rho =", opts_init.chem_rho
print "dt =", opts_init.dt

print "chem_switch = ", opts_init.chem_switch
print "coal_switch = ", opts_init.coal_switch
print "sedi_switch = ", opts_init.sedi_switch
print "src_switch = ", opts_init.src_switch

print "sstp_cond =", opts_init.sstp_cond
print "sstp_coal =", opts_init.sstp_coal
print "sstp_chem =", opts_init.sstp_chem 

print "n_sd_max =", opts_init.n_sd_max
print lgrngn.backend_t.OpenMP
print lgrngn.backend_t.CUDA
backend = lgrngn.backend_t.serial

opts = lgrngn.opts_t()
opts.sedi=0
opts.coal=0
opts.cond=0
opts.adve=1
print "adve =", opts.adve
print "sedi =", opts.sedi
print "cond =", opts.cond
print "coal =", opts.coal
print "src =", opts.src
print "chem_dsl =", opts.chem_dsl
print "chem_dcs =", opts.chem_dsc
print "chem_rct =", opts.chem_rct
print "RH_max =", opts.RH_max

opts.chem_gas = {
  lgrngn.chem_species_t.SO2  : 44,
  lgrngn.chem_species_t.O3   : 44,
  lgrngn.chem_species_t.H2O2 : 44
}
print "chem_gas[SO2] = ", opts.chem_gas[lgrngn.chem_species_t.SO2]
print "chem_gas = ", opts.chem_gas

# 0D (parcel)
print "0D"
rhod = arr_t([  1.])
th   = arr_t([300.])
rv   = arr_t([  0.01])

prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)
try: 
  prtcls.init(th, rv, rhod)
  raise Exception("multiple init call not reported!")
except:
  pass
prtcls.step_sync(opts, th, rv, rhod)
try:
  prtcls.step_sync(opts, th, rv, rhod)
  raise Exception("sync/async order mismatch not reported!")
except:
  pass
rain = prtcls.step_async(opts)
prtcls.step_sync(opts, th, rv)
prtcls.diag_dry_rng(0.,1.)
prtcls.diag_wet_rng(0.,1.)
prtcls.diag_dry_mom(1)
prtcls.diag_wet_mom(1)
prtcls.diag_all()
#prtcls.diag_chem(lgrngn.chem_species_t.OH)
prtcls.diag_sd_conc()
assert frombuffer(prtcls.outbuf()) == opts_init.sd_conc # parcel set-up

# 1D (periodic horizontal domain)
print "1D"
rhod = arr_t([  1.,   1.,   1.])
th   = arr_t([300., 300., 300.])
rv   = arr_t([   .01,  .01,  .01])
C    = arr_t([   .5,   .5,   .5,  .5])

opts_init.nx = 3 
opts_init.dx = 10
opts_init.x1 = opts_init.nx * opts_init.dx


print rhod.shape
print th.shape
print rv.shape
print C.shape

for it in range(2):
  prtcls = lgrngn.factory(backend, opts_init)
  if(it==0):
    prtcls.init(th, rv, rhod)
  else:  
    prtcls.init(th, rv, rhod, C)

  prtcls.step_sync(opts, th, rv)
  prtcls.step_async(opts)
  prtcls.step_sync(opts, th, rv, rhod)
  prtcls.step_async(opts)
  prtcls.step_sync(opts, th, rv, rhod, C)

  prtcls.diag_sd_conc()
  assert len(frombuffer(prtcls.outbuf())) == opts_init.nx
  print frombuffer(prtcls.outbuf())
  assert (frombuffer(prtcls.outbuf()) > 0).all()
  assert sum(frombuffer(prtcls.outbuf())) == opts_init.nx * opts_init.sd_conc

# 2D (periodic horizontal domain)
print "2D"
rhod = arr_t([[  1.,    1.   ],     [   1.,     1.  ]])
th   = arr_t([[300.,  300.   ],     [ 300.,   300.  ]])
rv   = arr_t([[   .01,   .01 ],     [    .01,    .01]])
Cx   = arr_t([[   .5,    .5], [    .5,     .5], [    .5,     .5]])
Cz   = arr_t([[   .0,    .0,   .0], [   .0,      .0,  .0]])

opts_init.nz = 2
opts_init.nx = 2
opts_init.dz = 10
opts_init.dx = 10
opts_init.z1 = opts_init.nz * opts_init.dz
opts_init.x1 = opts_init.nx * opts_init.dx

prtcls = lgrngn.factory(backend, opts_init)
try:
  prtcls.init(th, rv, rhod, Cx=Cx, Cy=Cz) 
  raise Exception("2D using Cx, Cy not reported!")
except:
  pass

for it in range(2):
  prtcls = lgrngn.factory(backend, opts_init)
  if(it==0):
    prtcls.init(th, rv, rhod)
  else:  
    prtcls.init(th, rv, rhod, Cx=Cx, Cz=Cz)

  prtcls.step_sync(opts, th, rv)
  prtcls.step_async(opts)
  prtcls.step_sync(opts, th, rv, rhod)
  prtcls.step_async(opts)
  prtcls.step_sync(opts, th, rv, rhod, Cx=Cx, Cz=Cz)
  prtcls.step_async(opts)
  prtcls.step_sync(opts, th, rv, rhod, Cx=Cx, Cz=Cz)

  prtcls.diag_sd_conc()
  assert len(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.nx
  print frombuffer(prtcls.outbuf())
  assert (frombuffer(prtcls.outbuf()) > 0).all()
  assert sum(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.nx * opts_init.sd_conc
  assert opts_init.nx == prtcls.opts_init.nx

#TODO: test profile vs. 2D array

# 3D
print "3D"
rhod = arr_t([rhod, rhod])
th   = arr_t([th,   th  ])
rv   = arr_t([rv,   rv  ])

opts_init.ny = 2
opts_init.dy = 10
opts_init.y1 = opts_init.ny * opts_init.dy

Cx = 0.5 * ones((opts_init.nx+1, opts_init.ny+0, opts_init.nz+0), dtype=float64) #TODO: these dimensions are not checked...
Cy = 0.5 * ones((opts_init.nx+0, opts_init.ny+1, opts_init.nz+0), dtype=float64)
Cz = zeros((opts_init.nx+0, opts_init.ny+0, opts_init.nz+1), dtype=float64)

for it in range(2):
  prtcls = lgrngn.factory(backend, opts_init)
  if(it==0):
    prtcls.init(th, rv, rhod)
  else:  
    prtcls.init(th, rv, rhod, Cx=Cx, Cy=Cy, Cz=Cz)

  prtcls.step_sync(opts, th, rv)
  prtcls.step_async(opts)
  prtcls.step_sync(opts, th, rv, rhod)
  prtcls.step_async(opts)
  prtcls.step_sync(opts, th, rv, rhod, Cx=Cx, Cy=Cy, Cz=Cz)
  prtcls.diag_sd_conc()

  assert len(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.nx * opts_init.ny
  print frombuffer(prtcls.outbuf())
  assert (frombuffer(prtcls.outbuf()) > 0).all()
  assert sum(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.nx * opts_init.ny * opts_init.sd_conc

#TODO: test profile vs. 3D array
