import sys
#try:
#  import boost.mpi
#except:
#  pass
sys.path.insert(0, "../../bindings/python/")
sys.path.insert(0, "../../../build/bindings/python/")

from libcloudphxx import lgrngn

from numpy import array as arr_t, frombuffer, repeat, zeros, float64, ones, isclose, asarray

from math import exp, log, sqrt, pi

n_tot  = 60e6
def lognormal(lnr):
  mean_r = .04e-6 / 2
  stdev  = 1.4
  return n_tot * exp(
    -pow((lnr - log(mean_r)), 2) / 2 / pow(log(stdev),2)
  ) / log(stdev) / sqrt(2*pi);

opts_init = lgrngn.opts_init_t()
kappa1 = .61
kappa2 = 1.28
kappa3 = 0.8
rd_insol = 0.5e-6
rho_stp = 1.2248
opts_init.dry_distros = {(kappa1, rd_insol): lognormal}
opts_init.kernel = lgrngn.kernel_t.geometric
opts_init.terminal_velocity = lgrngn.vt_t.beard76
opts_init.adve_scheme = lgrngn.as_t.euler
opts_init.dt = 1
opts_init.sd_conc = 64
opts_init.n_sd_max = int(1e6) # some space for tail SDs
opts_init.rng_seed = 396
opts_init.rng_seed_init = 456
opts_init.rng_seed_init_switch = True
opts_init.src_dry_distros = {kappa1:lognormal}
opts_init.rlx_dry_distros = {kappa1:[lognormal, [0,2], [0,opts_init.dz]]}
opts_init.src_sd_conc = 64
opts_init.src_z1 = opts_init.dz
opts_init.diag_incloud_time = True

print("nx = ", opts_init.nx)
print("ny = ", opts_init.ny)
print("nz = ", opts_init.nz)
print("dx = ", opts_init.dx)
print("dy = ", opts_init.dy)
print("dz = ", opts_init.dz)
print("x0 = ", opts_init.x0)
print("y0 = ", opts_init.y0)
print("z0 = ", opts_init.z0)
print("x1 = ", opts_init.x1)
print("y1 = ", opts_init.y1)
print("z1 = ", opts_init.z1)

print("sd_conc = ", opts_init.sd_conc)
print("RH_max = ", opts_init.RH_max)
print("rng_seed = ", opts_init.rng_seed)
print("rng_seed_init = ", opts_init.rng_seed_init)
print("kernel =", opts_init.kernel) 
print("sd_conc =", opts_init.sd_conc)
print("terminal_velocity =", opts_init.terminal_velocity)
print("adve_scheme =", opts_init.adve_scheme)
print("chem_rho =", opts_init.chem_rho)
print("dt =", opts_init.dt)

print("chem_switch = ", opts_init.chem_switch)
print("coal_switch = ", opts_init.coal_switch)
print("sedi_switch = ", opts_init.sedi_switch)
print("subs_switch = ", opts_init.subs_switch)
print("ice_switch = ", opts_init.ice_switch)
print("src_type = ", opts_init.src_type)

print("exact_sstp_cond = ", opts_init.exact_sstp_cond)
print("sstp_cond =", opts_init.sstp_cond)
print("sstp_coal =", opts_init.sstp_coal)
print("sstp_chem =", opts_init.sstp_chem) 

print("diag_incloud_time = ", opts_init.diag_incloud_time)
print("n_sd_max =", opts_init.n_sd_max)
print(lgrngn.backend_t.OpenMP)
print(lgrngn.backend_t.CUDA)
backend = lgrngn.backend_t.serial

opts = lgrngn.opts_t()
opts.sedi=0
opts.coal=0
opts.cond=0
opts.adve=1
print("adve =", opts.adve)
print("sedi =", opts.sedi)
print("subs =", opts.subs)
print("cond =", opts.cond)
print("coal =", opts.coal)
print("src =", opts.src)
print("chem_dsl =", opts.chem_dsl)
print("chem_dcs =", opts.chem_dsc)
print("chem_rct =", opts.chem_rct)
print("RH_max =", opts.RH_max)

opts.chem_gas = {
  lgrngn.chem_species_t.SO2  : 44,
  lgrngn.chem_species_t.O3   : 44,
  lgrngn.chem_species_t.H2O2 : 44
}
print("chem_gas[SO2] = ", opts.chem_gas[lgrngn.chem_species_t.SO2])
print("chem_gas = ", opts.chem_gas)

# --------- test runs -----------

# ----------
# 0D (parcel)
print("0D")
opts_init.sedi_switch = False
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
prtcls.step_async(opts)
prtcls.step_sync(opts, th, rv)
prtcls.diag_dry_rng(0.,1.)
prtcls.diag_wet_rng(0.,1.)
prtcls.diag_kappa_rng(0.,2.)
prtcls.diag_kappa_rng_cons(.5,1.5)
prtcls.diag_dry_mom(1)
prtcls.diag_wet_mom(1)
prtcls.diag_kappa_mom(1)
prtcls.diag_incloud_time_mom(1)
puddle = prtcls.diag_puddle()
print('puddle: ', puddle)
#prtcls.diag_chem(lgrngn.chem_species_t.OH)
prtcls.diag_all()
prtcls.diag_sd_conc()
print(frombuffer(prtcls.outbuf()))
assert frombuffer(prtcls.outbuf()) == opts_init.sd_conc # parcel set-up

# ----------
# 0D (parcel) with explicit calls to sync_in and step_cond 
print("0D")
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
try:
  prtcls.step_cond(opts, th, rv)
  raise Exception("sync_in/cond order mismatch not reported!")
except:
  pass
prtcls.sync_in(th, rv, rhod)
prtcls.step_cond(opts, th, rv)
prtcls.step_async(opts)
prtcls.step_sync(opts, th, rv)
prtcls.diag_dry_rng(0.,1.)
prtcls.diag_wet_rng(0.,1.)
prtcls.diag_dry_mom(1)
prtcls.diag_wet_mom(1)
prtcls.diag_kappa_mom(1)
prtcls.diag_incloud_time_mom(1)
puddle = prtcls.diag_puddle()
print('puddle: ', puddle)
#prtcls.diag_chem(lgrngn.chem_species_t.OH)
prtcls.diag_all()
prtcls.diag_sd_conc()
print(frombuffer(prtcls.outbuf()))
assert frombuffer(prtcls.outbuf()) == opts_init.sd_conc # parcel set-up


# ----------
# 0D with large_tail option
print("0D large tail")
opts_init.sd_conc_large_tail = 1
prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)

prtcls.step_sync(opts, th, rv)
prtcls.step_async(opts)
opts_init.sd_conc_large_tail = 0
prtcls.diag_all()
prtcls.diag_sd_conc()

assert len(frombuffer(prtcls.outbuf())) == 1
print(frombuffer(prtcls.outbuf()))
assert (frombuffer(prtcls.outbuf()) > 0).all()
assert sum(frombuffer(prtcls.outbuf())) >= opts_init.sd_conc



# ----------
# 0D const multi - number of SDs and number of particles
print("0D const multi")
sd_conc_old = opts_init.sd_conc
opts_init.sd_conc = 0
prtcls_per_cell = 2 * n_tot / rho_stp #rhod=1; 2* because of two distributions
opts_init.sd_const_multi = int(prtcls_per_cell / 64) 
prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)
prtcls.diag_all()
prtcls.diag_sd_conc()

prtcls.step_sync(opts, th, rv)
prtcls.step_async(opts)
prtcls.step_sync(opts, th, rv, rhod)
prtcls.step_async(opts)
prtcls.diag_all()
prtcls.diag_sd_conc()

assert len(frombuffer(prtcls.outbuf())) == 1
print(frombuffer(prtcls.outbuf()))
assert (frombuffer(prtcls.outbuf()) > 0).all()
sd_tot = frombuffer(prtcls.outbuf()).sum()
prtcls.diag_all()
prtcls.diag_wet_mom(0)
prtcls_tot = frombuffer(prtcls.outbuf()).sum()
assert ((prtcls_tot / sd_tot)  == opts_init.sd_const_multi)
opts_init.sd_const_multi = 0
opts_init.sd_conc = sd_conc_old



# ----------
# 0D dry_sizes init with two kappas
print("0D dry sizes")
opts_init.dry_distros = dict()
opts_init.dry_sizes = {(kappa1, rd_insol) : {1.e-6  : [30. * rho_stp, 15], 15.e-6 : [10. * rho_stp, 10]},
                       (kappa2, rd_insol) : {1.2e-6 : [20. * rho_stp, 10], 12.e-6 : [15. * rho_stp, 15]}}

sd_conc_old = opts_init.sd_conc
opts_init.sd_conc = 0
prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)

prtcls.diag_all()
prtcls.diag_sd_conc()
print(frombuffer(prtcls.outbuf()))
assert (frombuffer(prtcls.outbuf()) > 0).all()

sd_tot = frombuffer(prtcls.outbuf()).sum()
prtcls.diag_all()
prtcls.diag_wet_mom(0)
prtcls_tot = frombuffer(prtcls.outbuf()).sum()
print(frombuffer(prtcls.outbuf()))
assert (prtcls_tot == 75.) # 25 SDs have multiplicity = 2 and 25 have multiplicity = 1
assert (sd_tot == 50.) 

prtcls.diag_dry_rng(1e-6, 1.1e-6);
prtcls.diag_wet_mom(0)
n = frombuffer(prtcls.outbuf()).copy()
prtcls.diag_kappa_mom(1)
k = frombuffer(prtcls.outbuf())
print(n, k)
assert (n == 30 ).all()
assert isclose(k, n * kappa1, rtol=1e-20)

prtcls.diag_dry_rng(1.2e-6, 1.3e-6);
prtcls.diag_wet_mom(0)
n = frombuffer(prtcls.outbuf()).copy()
prtcls.diag_kappa_mom(1)
k = frombuffer(prtcls.outbuf())
print(n, k)
assert (n == 20 ).all()
assert isclose(k, n * kappa2, rtol=1e-20)

prtcls.diag_dry_rng(12e-6, 13e-6);
prtcls.diag_wet_mom(0)
n = frombuffer(prtcls.outbuf()).copy()
prtcls.diag_kappa_mom(1)
k = frombuffer(prtcls.outbuf())
print(n, k)
assert (n == 15 ).all()
assert isclose(k, n * kappa2, rtol=1e-20)

prtcls.diag_dry_rng(15e-6, 15.1e-6);
prtcls.diag_wet_mom(0)
n = frombuffer(prtcls.outbuf()).copy()
prtcls.diag_kappa_mom(1)
k = frombuffer(prtcls.outbuf())
print(n, k)
assert (n == 10 ).all()
assert isclose(k, n * kappa1, rtol=1e-20)

# go back to distros init
opts_init.sd_conc = sd_conc_old
opts_init.dry_sizes = dict()
opts_init.dry_distros = {(kappa1, rd_insol):lognormal, (kappa2, rd_insol):lognormal}



# ----------
# 0D dry_sizes + sd_conc init
print("0D dry_sizes + sd_conc")
opts_init.dry_sizes = {(kappa3, rd_insol) : {1.e-6 : [30. * rho_stp, 15], 15.e-6 : [10. * rho_stp, 5]}}

prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)

prtcls.diag_all()
prtcls.diag_sd_conc()
print(frombuffer(prtcls.outbuf()))
assert frombuffer(prtcls.outbuf())[0] == 84 # 64 from dry_distro and 20 from sizes

# go back to distros init
opts_init.dry_sizes = dict()



# ----------
# 0D dry_sizes + sd_conc + tail
print("0D dry_sizes + sd_conc + tail")
opts_init.dry_sizes = {(kappa3, rd_insol): {1.e-6 : [30. * rho_stp, 15], 15.e-6 : [10. * rho_stp, 5]}}
opts_init.sd_conc_large_tail = 1

prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)

prtcls.diag_all()
prtcls.diag_sd_conc()
print(frombuffer(prtcls.outbuf()))
assert frombuffer(prtcls.outbuf())[0] > 84 # 64 from dry_distro and 20 from sizes + tail

# go back to distros init
opts_init.sd_conc_large_tail = 0
opts_init.dry_sizes = dict()



# ----------
# 0D dry_sizes + const_multi init
print("0D dry_sizes + const_multi")
opts_init.dry_sizes = {(kappa3, rd_insol) : {1.e-6 : [30. * rho_stp, 15], 15.e-6 : [10. * rho_stp, 5]}}
opts_init.sd_conc = 0
prtcls_per_cell = 2 * n_tot / rho_stp #rhod=1; 2* because of two distributions
opts_init.sd_const_multi = int(prtcls_per_cell / 64) 

prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)

prtcls.diag_all()
prtcls.diag_sd_conc()
print(frombuffer(prtcls.outbuf()))
assert frombuffer(prtcls.outbuf())[0] == 84 # 64 from dry_distro and 20 from sizes

# go back to distros init
opts_init.sd_conc = sd_conc_old
opts_init.sd_const_multi = 0
opts_init.dry_sizes = dict()


# ----------
# 0D turb
print("0D turb")
eps = arr_t([  1e-4])
opts_init.turb_coal_switch=True
opts.turb_coal=True
prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)
prtcls.step_sync(opts, th, rv, diss_rate=eps)
prtcls.step_async(opts)

prtcls.diag_all()
prtcls.diag_sd_conc()
assert len(frombuffer(prtcls.outbuf())) == 1
print(frombuffer(prtcls.outbuf()))
assert (frombuffer(prtcls.outbuf()) > 0).all()
assert sum(frombuffer(prtcls.outbuf())) == 64

opts_init.turb_coal_switch=False
opts.turb_coal=False

# ----------
# 0D ice
print("0D ice")
th   = arr_t([263.])
opts_init.ice_switch = True
opts_init.coal_switch = False
opts.ice_nucl = True
prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)
prtcls.step_sync(opts, th, rv)
prtcls.step_async(opts)

prtcls.diag_all()
prtcls.diag_sd_conc()
assert len(frombuffer(prtcls.outbuf())) == 1
print(frombuffer(prtcls.outbuf()))
assert (frombuffer(prtcls.outbuf()) > 0).all()
assert sum(frombuffer(prtcls.outbuf())) == 64

opts_init.ice_switch=False
opts.ice_nucl=False
opts_init.coal_switch = True



# ----------
# 1D (periodic horizontal domain)
print("1D")
rhod = arr_t([  1.,   1.,   1.])
th   = arr_t([300., 300., 300.])
rv   = arr_t([   .01,  .01,  .01])
C    = arr_t([   .5, .5, .5, .5]) # 3 cells = 4 edges on which C needs to be defined

opts_init.nx = 3 
opts_init.dx = 10
opts_init.x1 = opts_init.nx * opts_init.dx


print(rhod.shape)
print(th.shape)
print(rv.shape)
print(C.shape)

for it in range(2):
  prtcls = lgrngn.factory(backend, opts_init)
  if(it==0):
    prtcls.init(th, rv, rhod)
  else:  
    prtcls.init(th, rv, rhod, Cx=C)

  prtcls.step_sync(opts, th, rv)
  prtcls.step_async(opts)
  prtcls.step_sync(opts, th, rv, rhod)
  prtcls.step_async(opts)
  prtcls.step_sync(opts, th, rv, rhod, C)

  prtcls.diag_all()
  prtcls.diag_sd_conc()
  assert len(frombuffer(prtcls.outbuf())) == opts_init.nx
  print(frombuffer(prtcls.outbuf()))
  assert (frombuffer(prtcls.outbuf()) > 0).all()
  assert sum(frombuffer(prtcls.outbuf())) == opts_init.nx * opts_init.sd_conc

# ----------
# 2D (periodic horizontal domain)
print("2D")
opts_init.sedi_switch = True
rhod = arr_t([[  1.,    1.   ],     [   1.,     1.  ]])
th   = arr_t([[300.,  300.   ],     [ 300.,   300.  ]])
rv   = arr_t([[   .01,   .01 ],     [    .01,    .01]])
Cx   = arr_t([[    .5,     .5], [    .5,     .5], [   .5,    .5]])
Cz   = arr_t([[   .0,    .0,   .0], [   .0,      .0,  .0]])

opts_init.nz = 2
opts_init.nx = 2
opts_init.dz = 10
opts_init.dx = 10
opts_init.z1 = opts_init.nz * opts_init.dz
opts_init.x1 = opts_init.nx * opts_init.dx

opts_init.w_LS = zeros(opts_init.nz)
opts_init.SGS_mix_len = ones(opts_init.nz)
opts_init.subs_switch = True

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

  prtcls.diag_all()
  prtcls.diag_sd_conc()
  assert len(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.nx
  print(frombuffer(prtcls.outbuf()))
  assert (frombuffer(prtcls.outbuf()) > 0).all()
  assert sum(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.nx * opts_init.sd_conc
  assert opts_init.nx == prtcls.opts_init.nx

#TODO: test profile vs. 2D array

print("2D turb")
eps   = arr_t([[   1e-4,  1e-4], [   1e-4,  1e-4]])
opts_init.turb_adve_switch=True
opts.turb_adve=True
opts_init.turb_cond_switch=True
opts.turb_cond=True
prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)
prtcls.step_sync(opts, th, rv, diss_rate=eps)
prtcls.step_async(opts)

prtcls.diag_all()
prtcls.diag_sd_conc()
assert len(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.nx
print(frombuffer(prtcls.outbuf()))
assert (frombuffer(prtcls.outbuf()) > 0).all()
assert sum(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.nx * opts_init.sd_conc
assert opts_init.nx == prtcls.opts_init.nx

opts_init.turb_adve_switch=False
opts.turb_adve=False
opts_init.turb_cond_switch=False
opts.turb_cond=False


# ----------
# 3D
print("3D")
rhod = arr_t([rhod, rhod])
th   = arr_t([th,   th  ])
rv   = arr_t([rv,   rv  ])

opts_init.ny = 2
opts_init.dy = 10
opts_init.y1 = opts_init.ny * opts_init.dy

Cx = 0.5 * ones((opts_init.nx+1, opts_init.ny+0, opts_init.nz+0), dtype=float64) #TODO: these dimensions are not checked...
Cy = 0.5 * ones((opts_init.nx,   opts_init.ny+1, opts_init.nz+0), dtype=float64)
Cz =      zeros((opts_init.nx,   opts_init.ny+0, opts_init.nz+1), dtype=float64)

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
  prtcls.diag_all()
  prtcls.diag_sd_conc()

  assert len(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.nx * opts_init.ny
  print(frombuffer(prtcls.outbuf()))
  assert (frombuffer(prtcls.outbuf()) > 0).all()
  assert sum(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.nx * opts_init.ny * opts_init.sd_conc

print("3D turb")
eps   = arr_t([eps, eps])
opts_init.turb_adve_switch=True
opts.turb_adve=False
opts_init.turb_cond_switch=True
opts.turb_cond=True
prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)
prtcls.step_sync(opts, th, rv, diss_rate=eps)
prtcls.step_async(opts)

prtcls.diag_all()
prtcls.diag_sd_conc()
assert len(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.nx * opts_init.ny
print(frombuffer(prtcls.outbuf()))
assert (frombuffer(prtcls.outbuf()) > 0).all()
assert sum(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.nx * opts_init.ny * opts_init.sd_conc

opts_init.turb_adve_switch=False
opts.turb_adve=False
opts_init.turb_cond_switch=False
opts.turb_cond=False

# 3D with large_tail option
print("3D large tail")
opts_init.sd_conc_large_tail = 1
prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)

prtcls.step_sync(opts, th, rv)
prtcls.step_async(opts)
opts_init.sd_conc_large_tail = 0
prtcls.diag_all()
prtcls.diag_sd_conc()

assert len(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.nx * opts_init.ny
print(frombuffer(prtcls.outbuf()))
assert (frombuffer(prtcls.outbuf()) > 0).all()
assert sum(frombuffer(prtcls.outbuf())) >= opts_init.nz * opts_init.nx * opts_init.ny * opts_init.sd_conc



# ----------
# 3D const multi - number of SDs and number of particles
print("3D const multi")
opts_init.sd_conc = 0
cell_vol = opts_init.dx * opts_init.dy * opts_init.dz
prtcls_per_cell = 2 * n_tot * cell_vol / rho_stp #rhod=1; 2* because of two distributions
opts_init.sd_const_multi = int(prtcls_per_cell / 64) 
n_cell = opts_init.nz * opts_init.nx * opts_init.ny
prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)
prtcls.diag_all()
prtcls.diag_sd_conc()

prtcls.step_sync(opts, th, rv)
prtcls.step_async(opts)
prtcls.step_sync(opts, th, rv, rhod)
prtcls.step_async(opts)
prtcls.step_sync(opts, th, rv, rhod, Cx=Cx, Cy=Cy, Cz=Cz)
prtcls.diag_all()
prtcls.diag_sd_conc()

assert len(frombuffer(prtcls.outbuf())) == n_cell
print(frombuffer(prtcls.outbuf()))
assert (frombuffer(prtcls.outbuf()) > 0).all()
sd_tot = frombuffer(prtcls.outbuf()).sum()
prtcls.diag_all()
prtcls.diag_wet_mom(0)
prtcls_tot = frombuffer(prtcls.outbuf()).sum()
assert ((prtcls_tot / sd_tot) * cell_vol  == opts_init.sd_const_multi)



# ----------
# 3D dry_sizes init
print("3D dry sizes")
opts_init.dry_distros = dict()
opts_init.dry_sizes = {(kappa1, rd_insol) : {1.e-6 : [30./ cell_vol * rho_stp, 30], 15.e-6 : [10. / cell_vol * rho_stp, 10]}}

prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)

prtcls.diag_all()
prtcls.diag_sd_conc()
print(frombuffer(prtcls.outbuf()))
assert (frombuffer(prtcls.outbuf()) > 0).all()

sd_tot = frombuffer(prtcls.outbuf()).sum()
prtcls.diag_all()
prtcls.diag_wet_mom(0)
prtcls_tot = frombuffer(prtcls.outbuf()).sum()
print(frombuffer(prtcls.outbuf()))
assert ((prtcls_tot / sd_tot) * cell_vol  == 1)

prtcls.diag_dry_rng(1e-6, 1.1e-6);
prtcls.diag_wet_mom(0)
print(frombuffer(prtcls.outbuf()))
assert (frombuffer(prtcls.outbuf()) == 30 / cell_vol).all()

prtcls.diag_dry_rng(15e-6, 15.1e-6);
prtcls.diag_wet_mom(0)
print(frombuffer(prtcls.outbuf()))
assert (frombuffer(prtcls.outbuf()) == 10 / cell_vol).all()



# ----------
# 3D dry_sizes + sd_conc init
print("3D dry_sizes + sd_conc")
opts_init.dry_distros = {(kappa1, rd_insol):lognormal, (kappa2, rd_insol):lognormal}
opts_init.dry_sizes = {(kappa1, rd_insol) : {1.e-6 : [30./ cell_vol * rho_stp, 15], 15.e-6 : [10. / cell_vol * rho_stp,  5]}}
opts_init.sd_conc = sd_conc_old
opts_init.sd_const_multi = 0

prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)

prtcls.diag_all()
prtcls.diag_sd_conc()
print(frombuffer(prtcls.outbuf()))
assert (frombuffer(prtcls.outbuf()) == 84).all() # 64 from dry_distro and 20 from sizes

# test if get_attr work and if kappas are set correctly
kappa = asarray(prtcls.get_attr("kappa"))
# assert (kappa[:(32*opts_init.nx*opts_init.ny*opts_init.nz)] == kappa2).all()
# assert (kappa[(32*opts_init.nx*opts_init.ny*opts_init.nz):] == kappa1).all()
n = 32 * opts_init.nx * opts_init.ny * opts_init.nz
assert (kappa[:n] == kappa1).all()
assert (kappa[n:2*n] == kappa2).all()
assert (kappa[2*n:] == kappa1).all()


# ----------
# 3D dry_sizes + sd_conc + tail
print("3D dry_sizes + sd_conc + tail")
opts_init.sd_conc_large_tail = 1

prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)


prtcls.diag_all()
prtcls.diag_sd_conc()
print(frombuffer(prtcls.outbuf()))
assert (frombuffer(prtcls.outbuf())[0] > 64 + 20).all() # 64 from dry_distro and 20 from sizes + tail


# go back to distros init
opts_init.sd_conc_large_tail = 0
opts_init.dry_sizes = dict()



# ----------
# 3D dry_sizes + const_multi init
print("3D dry_sizes + const_multi")
opts_init.dry_sizes = {(kappa1, rd_insol) : {1.e-6 : [30./ cell_vol * rho_stp, 15], 15.e-6 : [10. / cell_vol * rho_stp, 5]}}
opts_init.sd_conc = 0
prtcls_per_cell = 2 * n_tot * cell_vol / rho_stp #rhod=1; 2* because of two distributions
opts_init.sd_const_multi = int(prtcls_per_cell / 64) 

prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)

prtcls.diag_all()
prtcls.diag_sd_conc()
print(frombuffer(prtcls.outbuf()))
assert (frombuffer(prtcls.outbuf())[0] == 84).all() # 64 from dry_distro and 20 from sizes

# go back to distros init
opts_init.sd_conc = sd_conc_old
opts_init.sd_const_multi = 0
opts_init.dry_sizes = dict()


# ----------
# 3D ice
print("3D ice")
th_arr   = arr_t([[263.,  263.   ],     [263.,  263.  ]])
th   = arr_t([th_arr,   th_arr  ])
opts_init.ice_switch = True
opts_init.coal_switch = False
opts.ice_nucl = True
prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)
prtcls.step_sync(opts, th, rv)
prtcls.step_async(opts)

prtcls.diag_all()
prtcls.diag_sd_conc()
assert len(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.nx * opts_init.ny
print(frombuffer(prtcls.outbuf()))
assert (frombuffer(prtcls.outbuf()) > 0).all()
assert sum(frombuffer(prtcls.outbuf())) == opts_init.nz * opts_init.nx * opts_init.ny * opts_init.sd_conc

opts_init.ice_switch = False
opts_init.coal_switch = True
opts.ice_nucl = False
