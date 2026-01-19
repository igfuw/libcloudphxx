import sys
sys.path.insert(0, "../../bindings/python/")

from libcloudphxx import lgrngn

from numpy import array as arr_t, frombuffer, repeat, zeros, float64, ones, isclose, mean

from math import exp, log, sqrt, pi

n_tot  = 60e6
def lognormal(lnr):
  mean_r = .04e-6 / 2
  stdev  = 1.4
  return n_tot * exp(
    -pow((lnr - log(mean_r)), 2) / 2 / pow(log(stdev),2)
  ) / log(stdev) / sqrt(2*pi);

def check_kappa_conc(prtcls, eps):
  prtcls.diag_kappa_rng(0.,1.)
  prtcls.diag_wet_mom(0)
  res_n = mean(frombuffer(prtcls.outbuf()))
  print(res_n * rho_stp)
  assert isclose(res_n * rho_stp, n_tot, atol=0., rtol=eps),\
    "initialized number of particles of type kappa1 differs from the distribution"
  
  prtcls.diag_kappa_rng(1.,2.)
  prtcls.diag_wet_mom(0)
  res_n = mean(frombuffer(prtcls.outbuf()))
  print(res_n * rho_stp)
  assert isclose(res_n * rho_stp, n_tot, atol=0., rtol=eps),\
    "initialized number of particles of type kappa2 differs from the distribution"

opts_init = lgrngn.opts_init_t()
kappa1 = .61
kappa2 = 1.28
rd_insol = 0.
rho_stp = 1.2248
opts_init.dry_distros = {(kappa1, rd_insol):lognormal, (kappa2, rd_insol):lognormal}
opts_init.kernel = lgrngn.kernel_t.geometric
opts_init.terminal_velocity = lgrngn.vt_t.beard76
opts_init.dt = 1
opts_init.sd_conc = 64
opts_init.n_sd_max = 512
opts_init.rng_seed = 396
opts_init.src_z1 = opts_init.dz
opts_init.sedi_switch = False

backend = lgrngn.backend_t.serial

opts = lgrngn.opts_t()

# 0D
rhod = arr_t([  1.])
th   = arr_t([300.])
rv   = arr_t([  0.01])

prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)

check_kappa_conc(prtcls, 2e-2)

# 3D
opts_init.ny = 2
opts_init.dy = 10
opts_init.y1 = opts_init.ny * opts_init.dy

opts_init.nx = 2
opts_init.dx = 10
opts_init.x1 = opts_init.nx * opts_init.dx

opts_init.nz = 2
opts_init.dz = 10
opts_init.z1 = opts_init.nz * opts_init.dz

rhod = 1. * ones((opts_init.nx, opts_init.ny, opts_init.nz), dtype=float64)
th = 300. * ones((opts_init.nx, opts_init.ny, opts_init.nz), dtype=float64)
rv = 0.01 * ones((opts_init.nx, opts_init.ny, opts_init.nz), dtype=float64)

prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)

check_kappa_conc(prtcls, 5e-3)

# 3D const multi - number of SDs and number of particles
opts_init.sd_conc = 0
cell_vol = opts_init.dx * opts_init.dy * opts_init.dz
prtcls_per_cell = 2 * n_tot * cell_vol / rho_stp #rhod=1
opts_init.sd_const_multi = int(prtcls_per_cell / 64) 
n_cell = opts_init.nz * opts_init.nx * opts_init.ny
opts_init.n_sd_max = int(n_cell * prtcls_per_cell / opts_init.sd_const_multi) # 2* because of two distributions
prtcls = lgrngn.factory(backend, opts_init)
prtcls.init(th, rv, rhod)

check_kappa_conc(prtcls, 5e-3)
