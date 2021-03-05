import sys
sys.path.insert(0, "../../bindings/python/")
#sys.path.insert(0, "/mnt/local/pdziekan/usr/local/lib/python3/dist-packages/")

from libcloudphxx import lgrngn

import numpy as np 
from math import exp, log, sqrt, pi
from time import time

def lognormal(lnr):
  mean_r = .04e-6 / 2
  stdev  = 1.4
  n_tot  = 60e6
  return n_tot * exp(
    -pow((lnr - log(mean_r)), 2) / 2 / pow(log(stdev),2)
  ) / log(stdev) / sqrt(2*pi);

Opts_init = lgrngn.opts_init_t()
kappa = .61
Opts_init.dry_distros = {kappa:lognormal}
Opts_init.coal_switch = False
Opts_init.sedi_switch = False
Opts_init.sgs_adve = lgrngn.sgs_adve_t.ST_periodic
Opts_init.ST_eps = 1e-4 # [m2/s3]
Opts_init.ST_Lmax = 100 # [m]
Opts_init.ST_Lmin = 1e-3 # [m]
Opts_init.Nmodes = 10
Opts_init.Nwaves_max = 10
Opts_init.init_pair_separation=0.1 #[m]

Opts_init.dt = 0.1

Opts_init.nz = 1
Opts_init.nx = 1
Opts_init.ny = 1
Opts_init.dz = 1e0
Opts_init.dx = 1e0
Opts_init.dy = 1e0
Opts_init.x0 = 0;
Opts_init.y0 = 0;
Opts_init.z0 = 0;
Opts_init.z1 = Opts_init.nz * Opts_init.dz
Opts_init.x1 = Opts_init.nx * Opts_init.dx
Opts_init.y1 = Opts_init.ny * Opts_init.dy
Opts_init.periodic_topbot_walls = 1

Opts_init.rng_seed = int(time())
Opts_init.sd_conc = 4
Opts_init.n_sd_max = Opts_init.sd_conc * (Opts_init.nx * Opts_init.nz * Opts_init.ny)

Backend = lgrngn.backend_t.CUDA
#Backend = lgrngn.backend_t.serial

Opts = lgrngn.opts_t()
Opts.adve = False
Opts.sgs_adve = True
Opts.sedi = False
Opts.cond = False
Opts.coal = False
Opts.chem = False
Opts.rcyc = False

Rhod       =   1. * np.ones((Opts_init.nx, Opts_init.ny, Opts_init.nz))
Th         = 300. * np.ones((Opts_init.nx, Opts_init.ny, Opts_init.nz))
Rv         = 0.01 * np.ones((Opts_init.nx, Opts_init.ny, Opts_init.nz))

prtcls = lgrngn.factory(Backend, Opts_init)
prtcls.init(Th, Rv, Rhod)

sep = prtcls.diag_pair_separation_mean()
print(sep)

for it in range(100):
  prtcls.step_sync(Opts, Th, Rv, Rhod)
  prtcls.step_async(Opts)
  sep = prtcls.diag_pair_separation_mean()
  print(sep)
