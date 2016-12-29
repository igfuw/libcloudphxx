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
opts_init.coal_switch = False
opts_init.sedi_switch = False

Opts = lgrngn.opts_t()
Opts.adve = True
Opts.sedi = False
Opts.cond = False
Opts.coal = False
Opts.chem = False

for adve_scheme in [lgrngn.as_t.euler, lgrngn.as_t.implicit, lgrngn.as_t.pred_corr]:
  opts_init.adve_scheme = adve_scheme

  try:
    prtcls = lgrngn.factory(lgrngn.backend_t.OpenMP, opts_init)
  except:
    prtcls = lgrngn.factory(lgrngn.backend_t.serial, opts_init)

  prtcls.init(th, rv, rhod)


  prtcls.step_sync(Opts,th,rv,rhod)
  prtcls.step_async(Opts)
