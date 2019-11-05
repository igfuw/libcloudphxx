import sys
sys.path.insert(0, "../../bindings/python/")

from libcloudphxx import lgrngn
from math import exp, log, sqrt, pi
import numpy as np

opts_init = lgrngn.opts_init_t()
opts_init.dt = 1

rhod = 1. * np.ones((1,))
th = 300. * np.ones((1,))
rv = 0.01 * np.ones((1,))

def lognormal(lnr):
  mean_r = .04e-6 / 2
  stdev = 1.4
  n_tot = 60e6
  return n_tot * exp(
    -pow((lnr - log(mean_r)), 2) / 2 / pow(log(stdev),2)
  ) / log(stdev) / sqrt(2*pi);

kappa = .61

opts_init.dry_distros = {kappa:lognormal}

opts_init.sd_conc = 50
opts_init.n_sd_max = 50

opts_init.kernel = lgrngn.kernel_t.geometric
opts_init.sedi_switch = False

Opts = lgrngn.opts_t()
Opts.adve = False
Opts.sedi = False
Opts.cond = False
Opts.coal = True
Opts.chem = False

for vt_eq in [lgrngn.vt_t.beard76, lgrngn.vt_t.beard77, lgrngn.vt_t.beard77fast, lgrngn.vt_t.khvorostyanov_spherical, lgrngn.vt_t.khvorostyanov_nonspherical]:
  opts_init.terminal_velocity = vt_eq

  try:
    prtcls = lgrngn.factory(lgrngn.backend_t.OpenMP, opts_init)
  except:
    prtcls = lgrngn.factory(lgrngn.backend_t.serial, opts_init)

  prtcls.init(th, rv, rhod)


  prtcls.step_sync(Opts,th,rv,rhod)
  prtcls.step_async(Opts)
