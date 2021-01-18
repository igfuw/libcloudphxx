import sys
sys.path.insert(0, "../../bindings/python/")

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
Opts_init.turb_adve_switch = True

Opts_init.dt = 1

Opts_init.nz = 6
Opts_init.nx = 6
Opts_init.dz = 1
Opts_init.dx = 1
Opts_init.z1 = Opts_init.nz * Opts_init.dz
Opts_init.x1 = Opts_init.nx * Opts_init.dx

Opts_init.rng_seed = int(time())
Opts_init.sd_conc = 100
Opts_init.n_sd_max = Opts_init.sd_conc * (Opts_init.nx * Opts_init.nz)

Backend = lgrngn.backend_t.serial

Opts = lgrngn.opts_t()
Opts.adve = False
Opts.turb_adve = True
Opts.sedi = False
Opts.cond = False
Opts.coal = False
Opts.chem = False
Opts.rcyc = False

Rhod       =   1. * np.ones((Opts_init.nx, Opts_init.nz))
Th         = 300. * np.ones((Opts_init.nx, Opts_init.nz))
Rv         = 0.01 * np.ones((Opts_init.nx, Opts_init.nz))
diss_rate  = 1e-4 * np.ones((Opts_init.nx, Opts_init.nz))

prtcls = lgrngn.factory(Backend, Opts_init)
prtcls.init(Th, Rv, Rhod)

prtcls.diag_all()
prtcls.diag_sd_conc()
tab_in = np.copy(np.frombuffer(prtcls.outbuf()).reshape(Opts_init.nx, Opts_init.nz))
print("at init \n", tab_in)

for it in range(100):
  prtcls.step_sync(Opts, Th, Rv, Rhod, diss_rate = diss_rate)
  prtcls.step_async(Opts)

prtcls.diag_all()
prtcls.diag_sd_conc()
tab_out = np.copy(np.frombuffer(prtcls.outbuf()).reshape(Opts_init.nx, Opts_init.nz))
print("after 100s \n", tab_out)

assert(np.array_equal(tab_in,tab_out) == False) # turbulence should have moved SDs around
