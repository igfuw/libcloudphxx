import sys
sys.path.insert(0, "../../bindings/python/")
sys.path.insert(0, "../../../build/bindings/python/")

from libcloudphxx import lgrngn

import numpy as np 
from math import exp, log, sqrt, pi
from time import time

def lognormal(lnr):
  mean_r = 100e-6
  stdev  = 1.4
  n_tot  = 1e6
  return n_tot * exp(
    -pow((lnr - log(mean_r)), 2) / 2 / pow(log(stdev),2)
  ) / log(stdev) / sqrt(2*pi);

Opts_init = lgrngn.opts_init_t()
kappa = .61
Opts_init.dry_distros = {kappa:lognormal}
Opts_init.coal_switch = False
Opts_init.sedi_switch = True
Opts_init.terminal_velocity = lgrngn.vt_t.beard76

Opts_init.dt = 1

Opts_init.nz = 1
Opts_init.nx = 2
Opts_init.dz = 1
Opts_init.dx = 1
Opts_init.z1 = Opts_init.nz * Opts_init.dz
Opts_init.x1 = Opts_init.nx * Opts_init.dx

Opts_init.rng_seed = int(time())
Opts_init.sd_conc = 10000
Opts_init.n_sd_max = Opts_init.sd_conc * (Opts_init.nx * Opts_init.nz)
Opts_init.dev_count = 2


Opts = lgrngn.opts_t()
Opts.adve = False
Opts.sedi = True
Opts.cond = False
Opts.coal = False
Opts.chem = False
Opts.rcyc = False

Rhod =   1. * np.ones((Opts_init.nx, Opts_init.nz))
Th   = 300. * np.ones((Opts_init.nx, Opts_init.nz))
Rv   = 0.01 * np.ones((Opts_init.nx, Opts_init.nz))

try:
  prtcls = lgrngn.factory(lgrngn.backend_t.multi_CUDA, Opts_init)
except:
  prtcls = lgrngn.factory(lgrngn.backend_t.serial, Opts_init)

prtcls.init(Th, Rv, Rhod)

for it in range(10):
  prtcls.step_sync(Opts, Th, Rv, Rhod)
  prtcls.step_async(Opts)

puddle = prtcls.diag_puddle()

prtcls.diag_all()
prtcls.diag_sd_conc()
tab_out = np.copy(np.frombuffer(prtcls.outbuf()).reshape(Opts_init.nx, Opts_init.nz))

assert(tab_out[0][0] == 0.)

puddle_expected_per_cell = {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: 0.0, 7: 0.0, 8: 7.087802417148837e-05, 9: 5.630090090571395e-06}

for a in puddle:
  print puddle[a], Opts_init.nx * puddle_expected_per_cell[a]
  assert np.isclose(puddle[a], Opts_init.nx * puddle_expected_per_cell[a], atol=0., rtol=1e-4)
